#!/usr/bin/env python

# coding: utf-8
# Copyright (c) MoGroup at UMD.
# Distributed under the terms of the MIT License.
from __future__ import print_function

"""
This script is used to analyze VASP AIMD result. It currently supports two options:

1. Analyze diffusivity from a series vasprun.xml (or vasprun.xml.gz) files at one
    temperature. It will output diffusivity and corresponding error analysis result.

2. Fitting Arrhenius relationship from diffusivities at different temperatures. It
    will output activation energy and corresponding error analysis result. It can
    also be used to predict diffusivity/conductivity at different temperature
    from the fitted Arrehenius relationship. The conversion factor from diffusivity
    to conductivity is (z*e)^2*c/(K*T), details can be found in paper
"""
import os
import argparse
from argparse import RawTextHelpFormatter
import csv
from prettytable import PrettyTable
from aimd.diffusion import DiffusivityAnalyzer, ErrorAnalysisFromDiffusivityAnalyzer, \
    ArreheniusAnalyzer, get_conversion_factor
from pymatgen.core.periodic_table import Specie
from pymatgen.io.vasp import Poscar

citing_info = """
-----------------------------------------
Thanks for using this script. 
Please consider citing our papers: https://github.com/mogroupumd/aimd#citing.
-----------------------------------------
"""


def Analyze_VASP_MD(args):
    """
    Analyze diffusivity from a series vasprun.xml (or vasprun.xml.gz) files at one
            temperature
    :param args: please check main function for details of args
    :return:
    """
    vasprun_dirs = []
    for i in range(args.runs_start, args.runs_end + 1):
        if os.path.exists(os.path.join(args.folder_feature + str(i), 'vasprun.xml.gz')):
            vasprun_dirs.append(os.path.join(args.folder_feature + str(i), 'vasprun.xml.gz'))
        elif os.path.exists(os.path.join(args.folder_feature + str(i), 'vasprun.xml')):
            vasprun_dirs.append(os.path.join(args.folder_feature + str(i), 'vasprun.xml'))
        else:
            raise Exception("No vasprun.xml or vasprun.xml.gz in folder {}".
                            format(args.folder_feature + str(i)))

    # In analyzing Arrhenius relationship, it is required to provide charged specie. To keep consistent, I also
    # require charged specie, even it is not necessary
    specie = Specie.from_string(args.specie)
    da = DiffusivityAnalyzer.from_files(vasprun_dirs, str(specie.element), step_skip=args.step_skip,
                                        ncores=args.ncores,
                                        time_intervals_number=args.time_intervals_number,
                                        spec_dict={'lower_bound': args.lower_bound_in_a_square \
                                                                  * args.site_distance \
                                                                  * args.site_distance,
                                                   'upper_bound': args.upper_bound,
                                                   'minimum_msd_diff': args.minimum_msd_diff_in_a_square \
                                                                       * args.site_distance \
                                                                       * args.site_distance,
                                                   }
                                        )
    ea = ErrorAnalysisFromDiffusivityAnalyzer(da, site_distance=args.site_distance)
    if da.diffusivity > 0:  # The linear fitting succeed
        summary_info = ea.get_summary_dict(oxidized_specie=args.specie)
    # if the msd profile of the MD doesn't fulfill the fitting requirements,
    # da.diffusivity is set  to be negative
    else:
        summary_info = {"diffusion result": "MSD calculated from MD doesn't fulfill the fitting requirement",
                        "max msd": max(da.msd),
                        'msd': da.msd,
                        'dt': da.dt,
                        'msd_component': da.msd_component
                        }
        print("Output msd-dt into {}K_msd-dt.csv".format(int(da.temperature)))
        args.msd_file = "{}K_msd-dt.csv".format(int(da.temperature))

    # output
    print("=" * 40)
    print("Used vasprun.xml files")
    print("Start run: {}, end run: {}".format(vasprun_dirs[0], vasprun_dirs[-1]))
    print("=" * 40)

    # If drift larger than 3 angstrom, raise warning
    if 'drift_maximum' in summary_info.keys():
        if max(summary_info['drift_maximum']) > 3:
            print("{} WARNING {}".format('='*20, '='*20))
            print("The entire cell has significant drift, please check the MD data")
            print("Maximum drift distance in xyz (Angstrom): {}".format(summary_info['drift_maximum']))
            print("{} WARNING {}".format('=' * 20, '=' * 20))
    # results table
    header_result = ("Parameter", "Value")
    result_table = PrettyTable(header_result)
    result_table.align["Parameter"] = "l"
    for k, v in summary_info.items():
        if k not in ['msd', 'dt', 'msd_component']:
            result_table.add_row([k, str(v)])
    result_table.add_row(['composition', str(da.structure.composition)])
    print("Results table: ")
    print("Diffusivity unit: cm^2/s, Conductivity unit: mS/cm, Drift unit: Angstrom")
    print(result_table.get_string(sortby='Parameter'))
    print(citing_info)
    # whether output msd
    if args.msd_file:
        print("Output msd-dt into file: {}".format(args.msd_file))
        with open(args.msd_file, 'w') as fp:
            w_csv = csv.writer(fp, delimiter=',')
            data = [summary_info['dt'], summary_info['msd'],
                    summary_info['msd_component'][0], summary_info['msd_component'][1],
                    summary_info['msd_component'][2]]
            w_csv.writerows([["dt (fs)", "msd (A^2)", "msd_component_0", "msd_component_1", "msd_component_2"]])
            w_csv.writerows(zip(*data))


def Analyze_arrhenius(args):
    """
    Fitting Arrhenius relationship from diffusivities at different temperatures.
    :param args:  please check main function for details of args
    :return:
    """
    aa = ArreheniusAnalyzer.from_csv(args.D_T_csv_file)
    # results table
    header_result = ("Parameter", "Value")
    result_table = PrettyTable(header_result)
    result_table.align["Parameter"] = "l"
    result_table.add_row(["Ea (eV) ", aa.Ea])
    result_table.add_row(["Ea error (eV)", aa.Ea_error])

    if args.temperature:
        result_table.add_row(["D at {} K (cm^2/s)".format(args.temperature),
                              aa.predict_diffusivity(args.temperature)[0]])
        result_table.add_row(["D range at {} K (cm^2/s)".format(args.temperature),
                              aa.predict_diffusivity(args.temperature)[1]])
    if args.poscar:
        if not args.specie:
            raise Exception("Please provide diffusion specie, with charged information, like Li+.")
        if not args.temperature:
            raise Exception("Please provide target temperature, like 300")
        poscar = Poscar.from_file(args.poscar)
        struct = poscar.structure
        result_table.add_row(["Convertion factor from D (cm^2/s) to S (mS/cm)",
                              get_conversion_factor(struct, args.specie, args.temperature)])
        result_table.add_row(["Conductivity at {} K (mS/cm)".format(args.temperature),
                              aa.predict_conductivity(args.temperature, struct, args.specie)[0]])
        result_table.add_row(["Conductivity range at {} K (mS/cm)".format(args.temperature),
                              aa.predict_conductivity(args.temperature, struct, args.specie)[1]])
    print("=" * 40)
    if aa.diffusivity_errors:
        print("Considering diffusivity errors as weight during analysis.")
    else:
        print("No diffusivity errors are provided, use ordinary fitting.")
    print("=" * 40)
    print("Results table: ")
    print(result_table.get_string())

    if args.verbose:
        verbose_table = PrettyTable(("Patameter", "Value"))
        verbose_table.add_row(["1000/T (1/K)", aa.x])
        verbose_table.add_row(["log10(D) (cm^2/s)", aa.y])
        if aa.y_error:
            verbose_table.add_row(["Error bar of log10(D) (cm^2/s)", aa.y_error])
        verbose_table.add_row(["Slope in arrhenius relationship", aa.slope])
        verbose_table.add_row(["Standard deviation of slope in arrhenius relationship", aa.slope_sigma])
        verbose_table.add_row(["Intercept in arrhenius relationship", aa.intercept])
        verbose_table.add_row(["Standard deviation of intercept in arrhenius relationship", aa.intercept_sigma])
        print("All details in arrhenius fitting")
        print(verbose_table.get_string())

    print(citing_info)
    if args.plot:
        plt = aa.get_arrhenius_plot()
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="""
        This script is used to analyze VASP AIMD result. It currently supports two options:

        1. Analyze diffusivity from a series vasprun.xml (or vasprun.xml.gz) files at one
            temperature. It will output diffusivity and corresponding error analysis result.

        2. Fitting Arrhenius relationship from diffusivities at different temperatures. It
            will output activation energy and corresponding error analysis result. It can
            also be used to predict diffusivity/conductivity at different temperature
            from the fitted Arrehenius relationship. The conversion factor from diffusivity
            to conductivity is (z*e)^2*c/(K*T), details can be found in paper ...
        """,
        formatter_class=RawTextHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_diffusivity = subparsers.add_parser(
        "diffusivity",
        help="Analyze diffusivity from a series vasprun.xml (or vasprun.xml.gz) files"
    )
    parser_diffusivity.add_argument(
        "specie", type=str,
        help="The diffusing species. "
             "Pleae provide species with charge, like Li+, not Li. Useful for calculation "
             "of conductivity from diffusivity."
    )
    parser_diffusivity.add_argument(
        "folder_feature", type=str,
        help="The shared folder feature of all MD result folders. "
             "For example, you want to analyze ./RUN_10/vasprun.xml.gz to ./RUN_100/vasprun.xml.gz. "
             "Then you should provide ./RUN_ or just RUN_. "
             "If your folders are run_000 to run_100, please change to format folder_feature + index_number, "
             "like run_0, ... run_100"
    )
    parser_diffusivity.add_argument(
        "runs_start", type=int,
        help="The start index number of the folders you want to analyze. "
    )
    parser_diffusivity.add_argument(
        "runs_end", type=int,
        help="The end index number of the folders you want to analyze. "
             "The end folder is included in the analysis. "
             "For example, if you want to analyze vasprun in folders [RUN_10, ..., RUN_20], "
             "you should provide runs_end as 20."
    )
    parser_diffusivity.add_argument(
        "site_distance", type=float,
        help="The distance between neighboring sites for the diffusion specie. "
             "If there are many different paths, you can provide an averaged value"
    )

    parser_diffusivity.add_argument(
        "-msd", "--msd_file", dest="msd_file", type=str, default=None,
        help="Option to save msd-dt file in csv format. "
             "If a file name is not provided, it will not output the msd-dt data"
    )
    parser_diffusivity.add_argument(
        "--step_skip", type=int, default=10,
        help="How many steps to skip in reading MD data, default is 10."
    )
    parser_diffusivity.add_argument(
        "-dt", "--time_intervals_number", type=int, default=1000,
        help="Spec of the time interval (dt) setting, the number of time_intervals. Default is 1000. "
    )
    parser_diffusivity.add_argument(
        "-l", "--lower_bound_in_a_square", type=float, default=0.5,
        help="The lower fitting bound in unit of a^2, a is the site distance. "
             "Default is 0.5. Typical lower bound is in unit of A^2"
    )
    parser_diffusivity.add_argument(
        "-u", "--upper_bound", type=float, default=0.5,
        help="The upper fitting bound, unit is total time. "
             "Default is 0.5 t_total"
    )
    parser_diffusivity.add_argument(
        "-min", "--minimum_msd_diff_in_a_square", type=float, default=0.5,
        help="The minimum msd range for fitting range, in unit of a^2. "
             "The fitting range should cover a range of MSD, default is 0.5. "
    )
    parser_diffusivity.add_argument(
        "-n", "--ncores", type=int, default=1,
        help="Number of cores used to do the analysis, default is 1. "
    )
    parser_diffusivity.set_defaults(
        func=Analyze_VASP_MD
    )

    parser_arrhenius = subparsers.add_parser(
        "arrhenius",
        help="Arrhenius relationship from diffusivities at different temperatures. "
             "It will read csv file containing diffusivities, temperatures and diffusivities error (optional). "
             "The header of csv MUST be D,T,D_error"
             "You can also predict diffusivity at different temperatures. "
             "More over, if you can also provide structure by POSCAR and the diffusion specie, "
             "you can obtain the corresponding conductivity."
    )  # Add this line to arrhenius -h output

    parser_arrhenius.add_argument(
        "D_T_csv_file", type=str,
        help="The csv file containing temperatures, diffusivities, and diffusivities error (error is optional). "
             "The csv header must be D,T, D_error, and the columns must follow that form."
    )
    parser_arrhenius.add_argument(
        "-T", "--temperature", type=float, default=None,
        help="new temperature to predict diffusivity/conductivity"
    )
    parser_arrhenius.add_argument(
        "-p", "--poscar", type=str, default=None,
        help="The POSCAR file of structure, any structure file is okay, "
             "as long as the volume and number of diffusion specie are same."
    )
    parser_arrhenius.add_argument(
        "-s", "--specie", type=str, default=None,
        help="The interested diffusion specie. "
             "It must charge information like Li+, not Li."
    )
    parser_arrhenius.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="Whether output details of fitting"
    )
    parser_arrhenius.add_argument(
        "-pl", "--plot", default=False, action='store_true',
        help="Whether plot the arrhenius plot."
    )
    parser_arrhenius.set_defaults(
        func=Analyze_arrhenius
    )

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
