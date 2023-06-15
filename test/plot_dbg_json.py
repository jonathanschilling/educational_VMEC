#!/usr/bin/env python3

import os
import sys
import json
import argparse

import numpy as np
import matplotlib.pyplot as plt

def plot_dbg_json(json_filename, quantities, surfaces):

    # put plots into a folder "plots" next to the given JSON file
    json_folder = os.path.normpath(os.path.dirname(json_filename))
    plots_folder = os.path.join(json_folder, "plots")
    os.makedirs(plots_folder, exist_ok=True)

    # assume that the given JSON file was created by educational_VMEC
    # --> first item of filename is context in which the corresponding file was created
    # --> last item of filename before ".json" is extension of corresponding VMEC run
    json_basename = os.path.basename(json_filename)

    parts_by_underscore = json_basename.split("_")
    dump_context = parts_by_underscore[0]

    first_dot_pos = json_basename.find(".")
    last_dot_pos = json_basename.rfind(".")
    if first_dot_pos < 0 or last_dot_pos < 0:
        raise RuntimeError("extension not found in JSON filename")
    extension = json_basename[first_dot_pos+1:last_dot_pos]

    # read all data from JSON file given as first command-line argument
    data = None
    with open(json_filename, "r") as f:
        data = json.load(f)

    # check that if the user specified a subset of quantities to plot,
    # all requested quantities are present in the given JSON file
    if quantities is not None:
        for quantity in quantities:
            if quantity not in data:
                print(f"ERROR: requested quantity '{quantity}' not found in JSON file '{json_filename}'.")
                sys.exit(-1)

    # iterate over all quantities in JSON file
    plt.figure(figsize=(10,5))
    num_keys = len(data.keys())
    for index, key in enumerate(data.keys()):

        # if a subset of quantities to plot is specified,
        # skip all that are not in the user-specified list to be plotted
        if quantities is not None and key not in quantities:
            continue

        value = np.array(data[key])

        num_dimensions = len(np.shape(value))
        if num_dimensions == 1:

            print("plot quantity %s (%d/%d) from context %s in VMEC run %s"%(
                key, index, num_keys,
                dump_context,
                extension
            ))

            plt.clf()
            plt.plot(value, ".-")
            plt.title(f"{extension} {key}")
            plt.tight_layout()

            fig_filename = f"{extension}_{dump_context}_{key}.pdf"
            plt.savefig(os.path.join(plots_folder, fig_filename))

        elif num_dimensions == 2:

            print("plot quantity %s (%d/%d) from context %s in VMEC run %s"%(
                key, index, num_keys,
                dump_context,
                extension
            ))

            plt.clf()
            plt.imshow(value[:, :].T, origin='upper', cmap='jet')
            plt.colorbar()
            plt.xlabel("toroidal grid point index")
            plt.ylabel("poloidal grid point index")
            plt.title(f"{extension} {key}")
            plt.tight_layout()

            fig_filename = f"{extension}_{dump_context}_{key}.pdf"
            plt.savefig(os.path.join(plots_folder, fig_filename))

        elif num_dimensions == 3:

            # hint: in forces JSON debug output, all arrays are shaped as follows: [ns, nzeta, ntheta3]
            ns = value.shape[0]

            # iterate over all flux surfaces
            for js in range(ns):

                # if a subset of surfaces to plot is specified,
                # skip all that are not in the user-specified list to be plotted
                if surfaces is not None and js not in surfaces:
                    continue

                print("plot quantity '%s' (%d/%d) at flux surface %d/%d from context '%s' in VMEC run '%s'"%(
                    key, index, num_keys,
                    js, ns,
                    dump_context,
                    extension
                ))

                plt.clf()
                plt.imshow(value[js, :, :].T, origin='upper', cmap='jet')
                plt.colorbar()
                plt.xlabel("toroidal grid point index")
                plt.ylabel("poloidal grid point index")
                plt.title(f"{extension} {key} js={js}")
                plt.tight_layout()

                fig_filename = f"{extension}_{dump_context}_{key}_{ns}_{js}.pdf"
                plt.savefig(os.path.join(plots_folder, fig_filename))
        else:
            print(f"skip {key} as it is not a 1D, 2D or 3D array (has {num_dimensions} dimensions)")
            continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot debugging output from JSON output of educational_VMEC"
    )
    parser.add_argument("filename", type=str, help="The JSON file to inspect.")
    parser.add_argument(
        "--quantities",
        type=str,
        nargs="+",
        help="The quantities from the JSON file to plot (default: all).",
    )
    parser.add_argument(
        "--surfaces",
        type=int,
        nargs="+",
        help="The flux surface indices (0 is axis, ns-1 is LCFS) to plot for 3D arrays (default: all).",
    )
    args = parser.parse_args()
    plot_dbg_json(json_filename=args.filename,
                  quantities=args.quantities,
                  surfaces=args.surfaces)