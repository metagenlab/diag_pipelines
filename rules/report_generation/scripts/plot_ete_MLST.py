import pandas
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, faces, NodeStyle,StackedBarFace
import matplotlib.colors as pcolors
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import rgb2hex

norm = mpl.colors.Normalize(vmin=float(snakemake.config['virulence_percentage_identity_cutoff']), vmax=100)
cmap_red = cm.OrRd
cmap_blue = cm.Blues
cmap_green = cm.Greens
m_red = cm.ScalarMappable(norm=norm, cmap=cmap_red)
m_blue = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
m_green = cm.ScalarMappable(norm=norm, cmap=cmap_green)

best_tree = snakemake.input["best_tree"]

leaf2mlst= pandas.read_csv(snakemake.input["mlst"],
                           delimiter='\t',
                           names=["leaf","species","mlst","1","2","3","4","5","6","7"]).set_index("leaf").to_dict()["mlst"]

def get_spaced_colors(n):

    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


def plot_tree_barplot(tree_file,
                      taxon2mlst,
                      header_list):

    '''

    display one or more barplot

    :param tree_file:
    :param taxon2value_list:
    :param exclude_outgroup:
    :param bw_scale:
    :param barplot2percentage: list of bool to indicates if the number are percentages and the range should be set to 0-100

    :return:
    '''

    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl

    mlst_list = list(set(taxon2mlst.values()))
    mlst2color = dict(zip(mlst_list, get_spaced_colors(len(mlst_list))))
    mlst2color['-'] = 'white'

    if isinstance(tree_file, Tree):
       t1 = tree_file
    else:
       t1 = Tree(tree_file)

    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)


    tss = TreeStyle()
    value=1
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False

    cmap = cm.YlGnBu#YlOrRd#OrRd


    scale_list = []
    max_value_list = []


    for i, lf in enumerate(t1.iter_leaves()):

        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:
            # header

            col_add=0



            #lf.add_face(n, column, position="aligned")
            n = TextFace('MLST')
            n.margin_top = 1
            n.margin_right = 2
            n.margin_left = 2
            n.margin_bottom = 1
            n.rotation = 90
            n.inner_background.color = "white"
            n.opacity = 1.
            n.hz_align = 2
            n.vt_align = 2

            tss.aligned_header.add_face(n, col_add+1)

        try:
            #if lf.name in leaf2mlst or int(lf.name) in leaf2mlst:
            n = TextFace(' %s ' % taxon2mlst[int(lf.name)])
            n.inner_background.color = 'white'
            m = TextFace('  ')
            m.inner_background.color = mlst2color[taxon2mlst[int(lf.name)]]
        except:
            n = TextFace(' na ')
            n.inner_background.color = "grey"
            m = TextFace('    ')
            m.inner_background.color = "white"

        n.opacity = 1.
        n.margin_top = 2
        n.margin_right = 2
        n.margin_left = 0
        n.margin_bottom = 2

        m.margin_top = 2
        m.margin_right = 0
        m.margin_left = 2
        m.margin_bottom = 2

        lf.add_face(m, 0, position="aligned")
        lf.add_face(n, 1, position="aligned")


        n = TextFace(lf.name, fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)

    for n in t1.traverse():
       nstyle = NodeStyle()
       if n.support < 1:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)

    return t1, tss

tree, style = plot_tree_barplot(best_tree,
                                leaf2mlst,
                                ["MLST"])

output_svg = snakemake.output[0]
tree.render(output_svg, dpi=1000, h=400, tree_style=style)