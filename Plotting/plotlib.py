import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import argparse
from cycler import cycler

#
# Constants
#

COLORMAP = 'viridis'
THESIS_WIDTH_PT = 398.33864
GOLDEN_RATIO = (5**.5 + 1) / 2

#
# Colors
#

def extract_cmap_colors(n, colormap=COLORMAP):
    cmap = plt.get_cmap(colormap)
    return cmap(np.linspace(0,.9,n))

def ccode_to_rgb(code):
    tmp = code.strip('#')
    if len(tmp) != 6:
        raise RuntimeError("Invalid color code")
    return int(tmp[0:2], 16) / 255, int(tmp[2:4], 16) / 255, int(tmp[4:6], 16) / 255

COLOR_PALETTE = [
    ccode_to_rgb('#4F0E5A'),
    ccode_to_rgb('#EE5622'),
    ccode_to_rgb('#ECA72C'),
    ccode_to_rgb('#222211'),
    ccode_to_rgb('#245C93'),
]

def set_color_cycler(colors):
    mpl.rc('axes', prop_cycle=(cycler(color=colors)))

#
# Plotting setup
#

def plot_setup(no_pgf=False, color_cnt=3, colors=COLOR_PALETTE, mode="thesis"):

    # reset to default parameters
    mpl.rcdefaults()

    set_color_cycler(colors)

    # LaTeX setup
    if not no_pgf:
        if mode == 'talk':
            main_font = "Fira Sans"
            sans_font = "Fira Sans"
            math_font = "Fira Math"
        else:
            sans_font = None
            main_font = "STIX Two Text"
            math_font = "STIX Two Math"

        font_params = [f"\\setmainfont{{{main_font}}}", f"\\setmathfont{{{math_font}}}"]
        if sans_font is not None:
            font_params.append(f"\\setsansfont{{{sans_font}}}")


        if mode == 'talk':
            cmyk = False
        else:
            cmyk = True

        # see https://matplotlib.org/3.1.1/gallery/userdemo/pgf_preamble_sgskip.html
        mpl.use("pgf")
        tex = {
                "font.family": "serif",
                "text.usetex": True,
                "pgf.texsystem": "lualatex",
                "pgf.rcfonts": False,
                "pgf.preamble": "\n".join([
                    # color
                    f"\\usepackage[{'cmyk' if cmyk else ''}]{{xcolor}}",

                    # font
                    r"\usepackage{fontspec}",
                    r"\usepackage{unicode-math}",
                    *font_params,

                    # additional packages
                    r"\usepackage{circledsteps}",
                    r"\usepackage{siunitx}",
                    r"\usepackage{xfrac}",
                    r"\usepackage{mhchem}",
                ]),
                }
        mpl.rcParams.update(tex)

    # fontsize parameters
    medium = {'thesis': 11, 'talk': 7}[mode]
    small = {'thesis': 9, 'talk': 6}[mode]

    mpl.rc('font', size=medium)
    mpl.rc('legend', fontsize=medium)
    mpl.rc('axes', labelsize=medium, titlesize=medium)
    mpl.rc('xtick', labelsize=small)
    mpl.rc('ytick', labelsize=small)

    # line parameters
    mpl.rc('lines', markersize=4.0)

    if mode=='talk':
        mpl.rc('figure', facecolor="None", edgecolor="None")

#
# Size calculation of plots
#
def set_size(width="thesis", fraction=1, aspect=GOLDEN_RATIO):

    if width == "thesis":
        width = THESIS_WIDTH_PT

    fig_width_pt = width * fraction
    inches_per_pt = 1 / 72.27

    fig_width_in = fig_width_pt * inches_per_pt
    fig_height_in = fig_width_in / aspect

    return (fig_width_in, fig_height_in)

# Enables the generation of ticks that are multiples of some number see https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
#

class MultipleTick:

    def __init__(self, denominator=2, number=np.pi, latex=r'\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    @staticmethod
    def formatter_proc(denominator, number, latex):
        def gcd(a, b):
            while b:
                a, b = b, a%b
            return a
        def _multiple_formatter(x, pos):
            den = denominator
            num = np.int(np.rint(den*x/number))
            com = gcd(num,den)
            (num,den) = (int(num/com),int(den/com))
            if den==1:
                if num==0:
                    return r'$0$'
                if num==1:
                    return r'$%s$'%latex
                elif num==-1:
                    return r'$-%s$'%latex
                else:
                    return r'$%s%s$'%(num,latex)
            else:
                if num==1:
                    return r'$\frac{%s}{%s}$'%(latex,den)
                elif num==-1:
                    return r'$\frac{-%s}{%s}$'%(latex,den)
                else:
                    return r'$\frac{%s%s}{%s}$'%(num,latex,den)
        return _multiple_formatter

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(self.formatter_proc(self.denominator, self.number, self.latex))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--palette", action="store_true")
    group.add_argument("--cmap", action="store_true")
    args = parser.parse_args()

    colors = None
    if args.palette:
        colors = COLOR_PALETTE
    elif args.cmap:
        colors = extract_cmap_colors(20, COLORMAP)

    if colors:
        with open('../TeX/preamble/colors.tex', 'w') as f:
            for cnt, color in enumerate(colors):
                f.write(f"\\definecolor{{cmap{cnt}}}{{rgb}}{{{color[0]:.5f},{color[1]:.5f},{color[2]:.5f}}}\n")
