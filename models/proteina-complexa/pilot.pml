python
from math import dist
from pathlib import Path
from pymol import cmd

complexes = sorted(path for path in Path("structures").glob("job_*/*.pdb")
                   if not path.stem.endswith("_binder"))
for number, path in enumerate(complexes, 1):
    name = f"design_{number}"
    cmd.load(str(path), name)
    binder = f"{name} and chain B"
    ligand = f"{name} and chain A"
    cmd.hide("everything", name)
    cmd.dss(binder)
    cmd.show("cartoon", binder)
    cmd.color("marine", f"{binder} and ss H")
    cmd.color("violet", f"{binder} and ss S")
    cmd.color("slate", f"{binder} and not (ss H+S)")
    cmd.show("sticks", ligand)
    cmd.color("gold", ligand)
    clashes = f"clashes_{number}"
    for left in cmd.get_model(f"{ligand} and not hydro").atom:
        for right in cmd.get_model(f"{binder} and not hydro").atom:
            if dist(left.coord, right.coord) < 2:
                cmd.distance(clashes, f"{name} and index {left.index}",
                             f"{name} and index {right.index}")
    cmd.set("dash_color", "red", clashes)
    cmd.hide("labels", clashes)
    cmd.group(f"complex_{number}", f"{name} {clashes}")
    cmd.disable("all")
    cmd.enable(f"complex_{number}")
    cmd.enable(name)
    cmd.enable(clashes)
    cmd.orient(name)
    cmd.zoom(name, 3)
    cmd.scene(f"Design {number}", "store")
cmd.bg_color("white")
cmd.scene("Design 1", "recall")
python end
