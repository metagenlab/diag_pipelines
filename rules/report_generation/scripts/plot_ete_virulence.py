import pandas
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, faces, NodeStyle
import matplotlib.colors as pcolors
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import rgb2hex

norm = mpl.colors.Normalize(vmin=80, vmax=100)
cmap_red = cm.OrRd
cmap_blue = cm.Blues
cmap_green = cm.Greens
m_red = cm.ScalarMappable(norm=norm, cmap=cmap_red)
m_blue = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
m_green = cm.ScalarMappable(norm=norm, cmap=cmap_green)

tab_files = [pandas.read_csv(name, delimiter='\t') for name in snakemake.input["nr_blast_results"]]
best_tree = snakemake.input["best_tree"]
leaf2spa_typing= pandas.read_csv(snakemake.input["spa_typing"],
                                 delimiter='\t',
                                 names=["SAMPLENAME","EGENOMICS_SPA_TYPE", "RIDOM_SPA_TYPE"]).set_index("SAMPLENAME").to_dict()["RIDOM_SPA_TYPE"]

leaf2mlst= pandas.read_csv(snakemake.input["mlst"],
                           delimiter='\t',
                           names=["leaf","spacies","mlst","1","2","3","4","5","6","7"]).set_index("leaf").to_dict()["mlst"]
vf_table = pandas.read_csv(snakemake.input["vf_table"], delimiter='\t')
vf_list = vf_table["gene"] + '_' + vf_table["uniprot_accession"]

leaf_id2protein_id2identity = {}
for blast_table, blast_file in zip(tab_files, snakemake.input["nr_blast_results"]):
    leaf_id = blast_file.split('/')[1]
    protein_id2identity = blast_table.set_index("virulence_factor_ID").to_dict()["percentage_identity"]
    leaf_id2protein_id2identity[leaf_id]= protein_id2identity


def plot_ete_tree(tree_file,
                  ordered_queries,
                  leaf_id2protein_id2identity,
                  leaf_id2description=False,
                  leaf_id2mlst=False,
                  leaf_id2spa=False,
                  show_identity_values=True):


    t1 = Tree(tree_file)
    tss = TreeStyle()
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    head = True
    column_add = 0
    for lf in t1.iter_leaves():
        lf.branch_vertical_margin = 0
        # add MLST
        if leaf_id2mlst:
            if head:
                column_add+=1
                n = TextFace(' MLST ')
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.rotation = 270
                n.vt_align = 2
                n.hz_align = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                tss.aligned_header.add_face(n, 0)

            if lf.name in leaf2mlst:
                n = TextFace(' %s ' % leaf2mlst[lf.name])
            else:
                n = TextFace(' na ')
            n.opacity = 1.
            n.margin_top = 2
            n.margin_right = 2
            n.margin_left = 2
            n.margin_bottom = 2
            n.inner_background.color = "grey"
            lf.add_face(n, 0, position="aligned")
        # add spa typing
        if leaf_id2spa:
            if head:
                column_add+=1
                n = TextFace(' spa ')
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.rotation = 270
                n.vt_align = 2
                n.hz_align = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                tss.aligned_header.add_face(n, column_add-1)
            if lf.name in leaf_id2spa:
                n = TextFace(' %s ' % leaf_id2spa[lf.name])
            else:
                n = TextFace(' na ')
            n.opacity = 1.
            n.margin_top = 2
            n.margin_right = 2
            n.margin_left = 2
            n.margin_bottom = 2
            n.inner_background.color = "grey"
            lf.add_face(n, column_add-1, position="aligned")
        # loop to add virulence gene hits
        for column, protein_id in enumerate(ordered_queries):
            # draw labels at the top of each column
            if head:
                if show_identity_values:
                    n = TextFace(' %s ' % str(protein_id))
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 2
                    n.rotation = 270
                    n.vt_align = 2
                    n.hz_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, column+column_add)
                else:
                    n = TextFace(' %s ' % str(protein_id), fsize=6)
                    n.margin_top = 0
                    n.margin_right = 0
                    n.margin_left = 0
                    n.margin_bottom = 0
                    n.rotation = 270
                    n.vt_align = 2
                    n.hz_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    # lf.add_face(n, col, position="aligned")
                    tss.aligned_header.add_face(n, column+column_add)
            # draw column content
            if lf.name not in leaf_id2protein_id2identity:
                n = TextFace(' %s ' % str('   -   '))
                n.opacity = 1.
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "grey"
                lf.add_face(n, column+column_add, position="aligned")
            else:
                if protein_id in leaf_id2protein_id2identity[lf.name]:
                    identity_value = float(leaf_id2protein_id2identity[lf.name][protein_id])
                    color = rgb2hex(m_blue.to_rgba(identity_value))


                    if show_identity_values:
                        # report identity values in coloured boxes
                        # adapt box size depending the digit width
                        if str(identity_value) == '100.00' or str(identity_value) == '100.0':
                            identity_value = '100'
                            n = TextFace(" %s  " % identity_value)
                        else:
                            n = TextFace("%.2f" % round(float(identity_value), 2))
                        # color text to white for dark cells
                        if float(identity_value) > 95:
                            n.fgcolor = "white"
                        n.opacity = 1.
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                        n.inner_background.color = color
                        lf.add_face(n, column+column_add, position="aligned")
                    else:
                        # draw coloured boxes without text
                        n = TextFace('  ')
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 0
                        # n.color = color
                        n.inner_background.color = color
                        lf.add_face(n, column+column_add, position="aligned")
                else:
                    n = TextFace('  %s  ' % str('  -  '))
                    n.opacity = 1.
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 2
                    n.inner_background.color = "white"
                    lf.add_face(n, column+column_add, position="aligned")

        # end of first leaf: turn off header
        head = False

    # draw boostrap supports
    for n in t1.traverse():
        nstyle = NodeStyle()
        if n.support < 0.9:
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 6
            n.set_style(nstyle)
        else:
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 0
            n.set_style(nstyle)

    return t1, tss

t1, tss = plot_ete_tree(best_tree,
              vf_list,
              leaf_id2protein_id2identity,
              leaf_id2spa=leaf2spa_typing,
              leaf_id2mlst=leaf2mlst)

output_svg = snakemake.output[0]
t1.render(output_svg, dpi=1000, h=400, tree_style=tss)