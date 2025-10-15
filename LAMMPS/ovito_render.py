#!/usr/bin/env python3
"""
OVITO Python script for rendering LAMMPS hard disk trajectories
Can be used with OVITO Pro or ovitos command-line tool
For OVITO Free: Load trajectories manually through GUI
"""

try:
    from ovito.io import import_file, export_file
    from ovito.modifiers import ColorCodingModifier, AssignColorModifier
    from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
    from ovito.pipeline import Pipeline
    import ovito
    OVITO_AVAILABLE = True
except ImportError:
    OVITO_AVAILABLE = False
    print("OVITO Python module not available.")
    print("This script requires OVITO Pro or 'ovitos' Python package.")
    print("\nFor visualization, use OVITO Free GUI:")
    print("  1. Download from https://www.ovito.org/")
    print("  2. Open OVITO")
    print("  3. File → Load File → Select .lammpstrj file")
    print("  4. Adjust visualization settings manually")
    import sys
    sys.exit(0)

import argparse
import numpy as np

def setup_visualization(pipeline):
    """Configure visualization pipeline"""

    # Color particles by type
    color_mod = ColorCodingModifier(
        property='Particle Type',
        gradient=ColorCodingModifier.Rainbow()
    )
    pipeline.modifiers.append(color_mod)

    # Override colors: Type 1 (particles) = Blue, Type 2 (wall) = Red
    # Note: In OVITO scripting, colors are RGB tuples (0-1 range)

    return pipeline

def render_snapshot(pipeline, output_file, frame=0, size=(1920, 1080)):
    """Render a single snapshot"""

    pipeline.compute(frame)

    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (0, 0, 50)  # Camera position (adjust for 2D view)
    vp.camera_dir = (0, 0, -1)   # Looking down -z axis

    vp.render_image(
        filename=output_file,
        size=size,
        renderer=TachyonRenderer()
    )

    print(f"Rendered snapshot: {output_file}")

def render_animation(pipeline, output_file, fps=30, size=(1920, 1080)):
    """Render full animation"""

    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (0, 0, 50)
    vp.camera_dir = (0, 0, -1)

    vp.render_anim(
        filename=output_file,
        size=size,
        fps=fps,
        renderer=TachyonRenderer()
    )

    print(f"Rendered animation: {output_file}")

def main():
    if not OVITO_AVAILABLE:
        return

    parser = argparse.ArgumentParser(
        description='Render LAMMPS trajectory with OVITO'
    )
    parser.add_argument('trajectory', help='LAMMPS trajectory file (.lammpstrj)')
    parser.add_argument('--output', default='render.png',
                        help='Output file (PNG for snapshot, MP4 for animation)')
    parser.add_argument('--animation', action='store_true',
                        help='Render full animation instead of snapshot')
    parser.add_argument('--frame', type=int, default=0,
                        help='Frame number for snapshot (default: 0)')
    parser.add_argument('--fps', type=int, default=30,
                        help='Frames per second for animation (default: 30)')
    parser.add_argument('--width', type=int, default=1920,
                        help='Output width in pixels (default: 1920)')
    parser.add_argument('--height', type=int, default=1080,
                        help='Output height in pixels (default: 1080)')

    args = parser.parse_args()

    # Load trajectory
    print(f"Loading trajectory: {args.trajectory}")
    pipeline = import_file(args.trajectory)

    # Setup visualization
    pipeline = setup_visualization(pipeline)

    # Render
    size = (args.width, args.height)

    if args.animation:
        render_animation(pipeline, args.output, fps=args.fps, size=size)
    else:
        render_snapshot(pipeline, args.output, frame=args.frame, size=size)

    print("Rendering complete!")

if __name__ == '__main__':
    main()
