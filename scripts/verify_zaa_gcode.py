#!/usr/bin/env python3
"""
Verify ZAA (Z Anti-Aliasing) G-code output.
Checks that Z varies within individual layers (non-planar extrusion).

Usage:
    python3 verify_zaa_gcode.py output_zaa.gcode [output_no_zaa.gcode]
"""

import re
import sys
from collections import defaultdict

def parse_gcode(filepath):
    """Parse G-code and extract per-layer Z variation data."""
    layers = []
    current_layer_z = None
    current_layer_moves = []
    z_values_in_layer = set()

    layer_pattern = re.compile(r';.*[Ll]ayer.*?(\d+)')
    g1_pattern = re.compile(r'G1\s')
    z_pattern = re.compile(r'Z([\d.]+)')
    e_pattern = re.compile(r'E([\d.]+)')

    with open(filepath) as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()

            # Detect layer change by Z-only moves or layer comments
            if line.startswith(';LAYER_CHANGE'):
                if current_layer_z is not None and z_values_in_layer:
                    layers.append({
                        'nominal_z': current_layer_z,
                        'z_values': sorted(z_values_in_layer),
                        'num_moves': len(current_layer_moves),
                    })
                z_values_in_layer = set()
                current_layer_moves = []
                continue

            if line.startswith(';Z_HEIGHT:'):
                try:
                    current_layer_z = float(line.split(':')[1])
                except (ValueError, IndexError):
                    pass
                continue

            # Track G1 moves with Z and E (extrusion moves with Z = contoured)
            if g1_pattern.match(line):
                z_match = z_pattern.search(line)
                e_match = e_pattern.search(line)
                if z_match and e_match:
                    z = float(z_match.group(1))
                    z_values_in_layer.add(round(z, 4))
                    current_layer_moves.append((line_no, z))

    # Don't forget the last layer
    if current_layer_z is not None and z_values_in_layer:
        layers.append({
            'nominal_z': current_layer_z,
            'z_values': sorted(z_values_in_layer),
            'num_moves': len(current_layer_moves),
        })

    return layers


def analyze_layers(layers):
    """Analyze layers for ZAA activity."""
    total_layers = len(layers)
    contoured_layers = 0
    max_variation = 0

    print(f"Total layers parsed: {total_layers}")
    print(f"{'='*60}")

    for i, layer in enumerate(layers):
        z_vals = layer['z_values']
        nominal = layer['nominal_z']

        if len(z_vals) <= 1:
            continue

        z_min = min(z_vals)
        z_max = max(z_vals)
        variation = z_max - z_min

        if variation > 0.001:  # More than 1 micron variation
            contoured_layers += 1
            max_variation = max(max_variation, variation)

            if contoured_layers <= 10:  # Show first 10 contoured layers
                print(f"Layer {i:4d} | nominal Z={nominal:.3f} | "
                      f"Z range: [{z_min:.4f}, {z_max:.4f}] | "
                      f"variation: {variation:.4f}mm | "
                      f"{len(z_vals)} unique Z values | "
                      f"{layer['num_moves']} extrusion moves")

    print(f"{'='*60}")
    print(f"Contoured layers: {contoured_layers}/{total_layers} "
          f"({100*contoured_layers/total_layers:.1f}%)" if total_layers > 0 else "No layers found")
    print(f"Max Z variation within a layer: {max_variation:.4f}mm")

    if contoured_layers > 0:
        print(f"\n*** ZAA IS WORKING *** — found {contoured_layers} layers with variable Z")
    else:
        print(f"\n*** ZAA NOT DETECTED *** — no layers have variable Z within them")
        print("  Possible causes:")
        print("  - zaa_enabled not set to 1 in config")
        print("  - Model has no top surfaces (flat top?)")
        print("  - ContourZ raycast found no mesh surface variation")

    return contoured_layers > 0


def compare_gcodes(zaa_file, no_zaa_file):
    """Compare ZAA vs non-ZAA G-code."""
    print(f"\n{'='*60}")
    print(f"COMPARISON: ZAA vs non-ZAA")
    print(f"{'='*60}")

    # Count G1 lines with Z+E (extrusion moves that include Z)
    def count_xyz_extrusions(filepath):
        count = 0
        with open(filepath) as f:
            for line in f:
                if re.match(r'G1\s.*Z[\d.].*E[\d.]', line) or re.match(r'G1\s.*E[\d.].*Z[\d.]', line):
                    count += 1
        return count

    zaa_count = count_xyz_extrusions(zaa_file)
    no_zaa_count = count_xyz_extrusions(no_zaa_file)

    print(f"G1 moves with both Z and E (XYZ extrusion):")
    print(f"  ZAA enabled:  {zaa_count}")
    print(f"  ZAA disabled: {no_zaa_count}")
    print(f"  Difference:   +{zaa_count - no_zaa_count}")

    if zaa_count > no_zaa_count + 10:
        print(f"\n*** CONFIRMED *** — ZAA adds {zaa_count - no_zaa_count} XYZ extrusion moves")
    else:
        print(f"\n*** WARNING *** — ZAA doesn't seem to add significant XYZ moves")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <zaa_gcode> [no_zaa_gcode]")
        sys.exit(1)

    print(f"Analyzing ZAA G-code: {sys.argv[1]}")
    print()

    layers = parse_gcode(sys.argv[1])
    zaa_ok = analyze_layers(layers)

    if len(sys.argv) >= 3:
        compare_gcodes(sys.argv[1], sys.argv[2])

    sys.exit(0 if zaa_ok else 1)
