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

blast_files = [pandas.read_csv(name, delimiter='\t') for name in snakemake.input["blast_results"]]
best_tree = snakemake.input["best_tree"]

ordered_samples = snakemake.params["samples"]

sample2n_VFs = {}
for n, sample in enumerate(ordered_samples):
    sample2n_VFs[sample] = [len(blast_files[n])]

leaf2mlst= pandas.read_csv(snakemake.input["mlst"],
                           delimiter='\t',
                           names=["leaf","species","mlst","1","2","3","4","5","6","7"]).set_index("leaf").to_dict()["mlst"]
nr_vf_list = []
sample2VF2identity={}
for blast_file in blast_files:
    # samples/{sample}/virulence/VFDB_results_blast.tsv
    sample = blast_file.split('/')[1]
    sample2VF2identity[]
    with open(blast_file, 'r') as f:
        for row in f:
            data = row.split('\t')
            if row[1] not in nr_vf_list:
                nr_vf_list.append(nr_vf_list)
            sample2VF2identity[row[1]] = row[2]



def get_spaced_colors(n):

    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


def plot_tree_barplot(tree_file,
                      taxon2value_list_barplot,
                      header_list,
                      taxon2set2value_heatmap=False,
                      header_list2=False,
                      column_scale=True,
                      general_max=False,
                      barplot2percentage=False,
                      taxon2mlst=False):

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

    if taxon2mlst:
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


    if column_scale and header_list2:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        for column in header_list2:
            values = taxon2set2value_heatmap[column].values()

            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m

    cmap = cm.YlGnBu#YlOrRd#OrRd

    values_lists = taxon2value_list_barplot.values()

    scale_list = []
    max_value_list = []

    for n, header in enumerate(header_list):
        #print 'scale', n, header
        data = [float(i[n]) for i in values_lists]

        if barplot2percentage is False:
            max_value = max(data)#3424182#
            min_value = min(data) #48.23
        else:
            if barplot2percentage[n] is True:
                max_value = 100
                min_value = 0
            else:
                max_value = max(data)#3424182#
                min_value = min(data) #48.23
        norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
        m1 = cm.ScalarMappable(norm=norm, cmap=cmap)
        scale_list.append(m1)
        if not general_max:
            max_value_list.append(float(max_value))
        else:
            max_value_list.append(general_max)

    for i, lf in enumerate(t1.iter_leaves()):

        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:

            col_add=0

            if taxon2mlst:
                header_list=['MLST']+header_list

            for col, header in enumerate(header_list):

                #lf.add_face(n, column, position="aligned")
                n = TextFace(' ')
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


                n = TextFace('%s' % header)
                n.margin_top = 1
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.rotation = 270
                n.inner_background.color = "white"
                n.opacity = 1.
                n.hz_align = 2
                n.vt_align = 1
                tss.aligned_header.add_face(n, col_add)
                col_add+=2

            if header_list2:
                for col, header in enumerate(header_list2):
                    n = TextFace('%s' % header)
                    n.margin_top = 1
                    n.margin_right = 20
                    n.margin_left = 2
                    n.margin_bottom = 1
                    n.rotation = 270
                    n.hz_align = 2
                    n.vt_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col+col_add)

        if taxon2mlst:

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
            col_add=2
        else:
            col_add=0

        try:
            val_list = taxon2value_list_barplot[lf.name]
        except:
            if not taxon2mlst:
                val_list = ['na'] * len(header_list)
            else:
                val_list = ['na'] * (len(header_list)-1)

        for col, value in enumerate(val_list):

            # show value itself
            try:
                n = TextFace('  %s  ' % str(value))
            except:
                n = TextFace('  %s  ' % str(value))
            n.margin_top = 1
            n.margin_right = 5
            n.margin_left = 10
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.

            lf.add_face(n, col_add, position="aligned")
            # show bar
            try:
                color = rgb2hex(scale_list[col].to_rgba(float(value)))
            except:
                color='white'
            try:
                percentage = (value/max_value_list[col])*100
                #percentage = value
            except:
                percentage = 0
            try:
                maximum_bar = ((max_value_list[col]-value)/max_value_list[col])*100
            except:
                maximum_bar = 0
            #maximum_bar = 100-percentage
            b = StackedBarFace([percentage,
                                maximum_bar],
                                width=100, height=10, colors=[color, "white"])
            b.rotation= 0
            b.inner_border.color = "grey"
            b.inner_border.width = 0
            b.margin_right = 15
            b.margin_left = 0
            lf.add_face(b, col_add+1, position="aligned")
            col_add+=2


        if taxon2set2value_heatmap:
            shift = col+col_add+1

            i = 0
            for col, col_name in enumerate(header_list2):
                try:
                    value = taxon2set2value_heatmap[col_name][lf.name]
                except:
                    try:
                        value = taxon2set2value_heatmap[col_name][int(lf.name)]
                    except:
                        value = 0

                if int(value) >0:
                    if int(value)>9:
                        n = TextFace(' %i ' % int(value))
                    else:
                        n = TextFace(' %i   ' % int(value))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.fgcolor = "white"
                    n.inner_background.color = rgb2hex(column2scale[col_name].to_rgba(float(value)))#"orange"
                    n.opacity = 1.
                    lf.add_face(n, col+col_add, position="aligned")
                    i+=1
                else:
                    n = TextFace('  ') #% str(value))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.inner_background.color = "white"
                    n.opacity = 1.

                    lf.add_face(n, col+col_add, position="aligned")


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
                                sample2n_VFs,
                                ["N hits"],
                                taxon2mlst=leaf2mlst)

output_svg = snakemake.output[0]
tree.render(output_svg, dpi=1000, h=400, tree_style=style)