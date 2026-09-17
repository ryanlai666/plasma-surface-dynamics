"""Render real kMC states as portable, GitHub-compatible GIFs.

Run: python -m plasma_surface.animate
The square layout is a display arrangement of independent surface columns.
"""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager
from .model import Parameters, ale_recipe, kmc
from .io import manifest, write_csv, write_json

BG = "#101b2e"
PANEL = "#18263d"
TEXT = "#edf3fb"
MUTED = "#a9bad2"
BARE = "#44c9c0"
MODIFIED = "#f1be53"
ETCHED = "#f074a5"
GRID = "#32475f"


def sample_trajectory(phases, p, *, cycles=3, side=16, seed=42, frames=121):
    if not isinstance(side, int) or side < 2 or not isinstance(frames, int) or frames < 2:
        raise ValueError("side and frames must be integers >= 2")
    events = []
    history = kmc(phases, p, cycles=cycles, sites=side*side, seed=seed, observer=events.append)
    times = np.array([event["time_s"] for event in events])
    requested = np.linspace(0, history[-1]["time_s"], frames)
    indices = np.searchsorted(times, requested, side="right")-1
    snapshots = [dict(events[int(i)], time_s=float(t)) for i,t in zip(indices,requested)]
    return history, snapshots


def _font(size):
    return ImageFont.truetype(font_manager.findfont("DejaVu Sans"), size)


def render_frame(name, snapshots, index, p, *, side, depth_layers, max_removal, phases, cycles):
    image = Image.new("RGB", (1000, 680), BG)
    draw = ImageDraw.Draw(image)
    snap = snapshots[index]
    modified = snap["modified"].reshape(side,side)
    heights = snap["heights"].reshape(side,side)
    previous = snapshots[max(0,index-1)]["heights"].reshape(side,side)
    just_etched = heights < previous
    title, body, small = _font(26), _font(16), _font(13)
    def text(x,y,value,fill=TEXT,font=body):
        draw.text((x,y),value,font=font,fill=fill)
    text(28,18,name+"  |  Surface kMC",font=title)
    text(28,55,"Illustrative rates | independent columns | no atom identities or chemical bonds",MUTED,small)
    text(738,22,f"Cycle {snap['cycle']} / {cycles}",MODIFIED)
    text(738,47,f"t = {snap['time_s']:5.2f} s   |   {snap['phase'].upper()}",TEXT,small)

    # Cycle timeline: proportional phase widths, with an absolute-time cursor.
    x0,x1,y0,y1=28,972,89,116
    duration=sum(phase.duration_s for phase in phases)
    total=duration*cycles
    colors={"modify":"#92722c","purge":"#304560","remove":"#27746e"}
    elapsed=0
    for cycle in range(cycles):
        for phase in phases:
            left=x0+(x1-x0)*elapsed/total
            elapsed+=phase.duration_s
            right=x0+(x1-x0)*elapsed/total
            draw.rectangle((left,y0,right-1,y1),fill=colors[phase.name])
            if right-left>65:
                text(left+6,y0+4,phase.name,TEXT,small)
    cursor=x0+(x1-x0)*snap['time_s']/total
    draw.line((cursor,y0-5,cursor,y1+5),fill=TEXT,width=3)

    draw.rounded_rectangle((28,139,424,485),radius=12,fill=PANEL)
    draw.rounded_rectangle((442,139,972,485),radius=12,fill=PANEL)
    text(48,152,"Top view: one tile = one column")
    text(463,152,"Cross-section: highlighted row")
    step=17
    gx,gy=83,195
    row=side//2
    for iy in range(side):
        for ix in range(side):
            left,top=gx+ix*step,gy+iy*step
            draw.rounded_rectangle((left,top,left+14,top+14),radius=2,
                                   fill=MODIFIED if modified[iy,ix] else BARE)
            if just_etched[iy,ix]:
                draw.rectangle((left-1,top-1,left+15,top+15),outline=ETCHED,width=2)
    draw.rectangle((gx-6,gy+row*step-4,gx+side*step+2,gy+(row+1)*step+2),outline=TEXT,width=2)
    text(75,466,"16 x 16 schematic site arrangement",MUTED,small)

    # Vertical coordinate is actual integer column height, with a fixed scale.
    sx,sy,bar_width=503,205,26
    layer_px=49
    bottom=sy+depth_layers*layer_px
    for depth in range(depth_layers+1):
        yy=sy+depth*layer_px
        draw.line((sx-6,yy,sx+side*bar_width,yy),fill=GRID,width=1)
        text(460,yy-8,f"{-depth}",MUTED,small)
    for ix in range(side):
        xx=sx+ix*bar_width
        top=sy-int(heights[row,ix])*layer_px
        draw.rectangle((xx,top,xx+bar_width-3,bottom),fill="#36516c")
        for yy in np.arange(top+layer_px,bottom,layer_px):
            draw.line((xx,yy,xx+bar_width-3,yy),fill=PANEL,width=2)
        draw.rectangle((xx,top,xx+bar_width-3,top+6),fill=MODIFIED if modified[row,ix] else BARE)
        if just_etched[row,ix]:
            draw.line((xx,top-4,xx+bar_width-3,top-4),fill=ETCHED,width=3)
    text(462,433,"Height in removal increments; zero = initial surface",MUTED,small)
    text(462,454,f"1 increment = {p.layer_nm:g} nm (illustrative conversion)",MUTED,small)

    for x,color,label in [(30,BARE,"Bare site"),(194,MODIFIED,"Modified site"),(395,ETCHED,"Etched since prior frame")]:
        draw.rectangle((x,502,x+13,515),fill=color)
        text(x+21,499,label,TEXT,small)
    coverage=float(snap['modified'].mean())
    removal=-float(snap['heights'].mean())*p.layer_nm
    text(712,496,f"Events: {snap['events']}",MUTED,small)
    text(712,516,f"Mean removal: {removal:.3f} nm",TEXT,small)

    def chart(box, values, ymax, label, color):
        l,t,r,b=box
        draw.line((l,t,l,b,r,b),fill=GRID,width=1)
        points=[]
        for frame,v in zip(snapshots[:index+1],values):
            point=(l+(r-l)*frame['time_s']/total,b-(b-t)*v/ymax)
            if points: points.append((point[0],points[-1][1]))
            points.append(point)
        if len(points)>1: draw.line(points,fill=color,width=2)
        if points:
            x,y=points[-1];draw.ellipse((x-3,y-3,x+3,y+3),fill=color)
        text(l,t-22,label,TEXT,small)
        text(l,b+5,"0 s",MUTED,small)
        text(r-37,b+5,f"{total:g} s",MUTED,small)
    coverages=[float(s['modified'].mean()) for s in snapshots[:index+1]]
    removals=[-float(s['heights'].mean())*p.layer_nm for s in snapshots[:index+1]]
    chart((40,568,420,630),coverages,1,f"Modified coverage: {coverage:.0%}   [scale 0-100%]",MODIFIED)
    chart((515,568,950,630),removals,max_removal,f"Mean removal   [scale 0-{max_removal:.3f} nm]",BARE)
    text(28,657,"Actual seeded Gillespie states; uniformly sampled simulation time. Visualization does not add lateral interactions.",MUTED,_font(12))
    return image


def run(output="docs/animations"):
    root=Path(__file__).resolve().parents[1]
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    cards=json.loads((root/"configs/materials.json").read_text(encoding="utf-8"))['cards']
    labels={"Si_reference":"Si", "SiN0.8_hypothesis":"SiN0.8", "SiN1.0_hypothesis":"SiN1.0", "Si3N4_hypothesis":"Si3N4"}
    phases=ale_recipe(energy_ev=35, dose_s=2, ion_s=4)
    cases=[]
    for card in cards:
        p=replace(Parameters(),**card['parameters'])
        history,snapshots=sample_trajectory(phases,p)
        # Exact parity confirms observing the trajectory did not alter random draws.
        assert history==kmc(phases,p,cycles=3,sites=256,seed=42)
        assert int(snapshots[-1]['heights'].sum())==snapshots[-1]['deposited']-snapshots[-1]['removed']
        cases.append((card,p,history,snapshots))
    depth=max(4,1-min(int(s['heights'].min()) for _,_,_,ss in cases for s in ss))
    ymax=max(.408,max(-float(s['heights'].mean())*p.layer_nm for _,p,_,ss in cases for s in ss)*1.1)
    summary=[]
    for card,p,history,snapshots in cases:
        label=labels[card['name']]
        print('Rendering',label,flush=True)
        frames=[render_frame(label,snapshots,i,p,side=16,depth_layers=depth,max_removal=ymax,phases=phases,cycles=3)
                for i in range(len(snapshots))]
        # One shared palette across every frame prevents color flicker.
        palette=frames[0].quantize(colors=128)
        gif_frames=[f.quantize(palette=palette,dither=Image.Dither.NONE) for f in frames]
        filename=label.lower().replace('.','p')+'.gif'
        gif_frames[0].save(output/filename,save_all=True,append_images=gif_frames[1:],duration=100,loop=0,disposal=2,optimize=False)
        frames[len(frames)//2].save(output/(label.lower().replace('.','p')+'_preview.png'))
        stats=[]
        for snap in snapshots:
            stats.append(dict(time_s=snap['time_s'],cycle=snap['cycle'],phase=snap['phase'],events=snap['events'],
                              coverage=float(snap['modified'].mean()),net_removed_nm=-float(snap['heights'].mean())*p.layer_nm))
        write_csv(output/(label.lower().replace('.','p')+'_frames.csv'),stats)
        with Image.open(output/filename) as gif:
            assert gif.n_frames==len(snapshots)
            durations=[]
            for i in range(gif.n_frames):
                gif.seek(i);durations.append(gif.info['duration'])
            assert sum(durations)==12100
        summary.append(dict(material=label,filename=filename,frames=len(snapshots),playback_duration_s=12.1,
                            simulation_duration_s=history[-1]['time_s'],final_net_removed_nm=history[-1]['net_removed_nm'],
                            events=history[-1]['events'],parameters=asdict(p),gif_sha256=hashlib.sha256((output/filename).read_bytes()).hexdigest()))
    write_json(output/'manifest.json',manifest(dict(backend='python',seed=42,sites=256,cycles=3,
        frames=121,frame_duration_ms=100,phase_schedule=[asdict(s) for s in phases],
        depth_scale_layers=depth,removal_scale_nm=ymax,scenarios=summary,
        interpretation='Independent-column display, not an atomistic Si/SiNx lattice. Nitride rates hypothetical.')))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    run()
