"""
Script to draw any non-canonical graph.
The main function just draws all the recurrent interior loops
(in Data/clusters_IL.pickle) in the folder Data/Figures
"""


import os
import pickle


import matplotlib
import matplotlib.patches as mpatches
import networkx as nx
from matplotlib import pyplot as plt

PATH_FIGURES = os.path.join('..', 'Data', 'Figures')
PATH_CLUSTERS = os.path.join('..' 'Data', 'clusters_IL.pickle')


def draw_graph(g0, title="", savefig=None):
    labels = {}
    elabels = {}

    pos = nx.circular_layout(g0)
    nodes = nx.draw_networkx_nodes(g0, pos, nodelist=g0.nodes(), node_color='lightgrey', node_size=500, linewidth=2,
                                   alpha=0.99)
    nodes.set_edgecolor("black")
    for index, i in enumerate(list(g0.nodes(data=True))):
        labels[i[0]] = i[1]["nt"]

    nx.draw_networkx_labels(g0, pos, labels, font_size=10)

    NC_edges = []
    C_edges = []
    backbone = []
    stacking = []
    for i in g0.edges(data=True):
        if i[0] > i[1]:
            continue
        i = i[0:2]
        label = g0.get_edge_data(*i)['label'].upper()
        if label not in ["B53", "S55", "S33"]:
            if label == "CWW":
                elabels[i] = r"●"
            if label == "TWW":
                elabels[i] = r"○"
            if label == "CSS":
                elabels[i] = r"▲"
            if label == "TSS":
                elabels[i] = r"∆"
            if label == "CHH":
                elabels[i] = r"■"
            if label == "THH":
                elabels[i] = r"□"

            if label == "CHW":
                elabels[i] = r"■●"
            if label == "THW":
                elabels[i] = r"□○"
            if label == "CWH":
                elabels[i] = r"●■"
            if label == "TWH":
                elabels[i] = r"○□"
            if label == "TSH":
                elabels[i] = r"∆□"
            if label == "CSH":
                elabels[i] = r"▲■"
            if label == "THS":
                elabels[i] = r"□∆"
            if label == "CHS":
                elabels[i] = r"■▲"
            if label == "CWS":
                elabels[i] = r"■●"
            if label == "TWS":
                elabels[i] = r"○∆"
            if label == "CSW":
                elabels[i] = r"▲●"
            if label == "TSW":
                elabels[i] = r"∆○"
        #elif label in ["S35", "S53", "S55", "S33"]:
        #    elabels[i] = r"s     "
        else:
            elabels[i] = r""
        if label.upper() == "CWW":
            C_edges.append(i)
        if label.upper() == "B53":
            backbone.append(i)
        if label.upper() in ["S55", "S33", "S35", "S53"]:
            stacking.append(i)
        if label.upper() not in ["S55", "S33", "S35", "S53", "B53", "CWW"]:
            NC_edges.append(i)
    #print(NC_edges)
    #print(elabels)
    nx.draw_networkx_edge_labels(g0, pos, elabels, font_size=10)
    nx.draw_networkx_edges(g0, pos, edgelist=NC_edges, edge_color='purple', width=4, arrows=False)
    nx.draw_networkx_edges(g0, pos, edgelist=C_edges, edge_color='green', width=4, arrows=False)
    nx.draw_networkx_edges(g0, pos, edgelist=backbone, edge_color='black', width=3, arrowsize=25)
    nx.draw_networkx_edges(g0, pos, edgelist=stacking, edge_color='orange', width=1, arrows=False)
    plt.title(title, fontsize=24)
    #print("done")
    plt.axis("off")
    NCP = mpatches.Patch(color="purple", label="Non-canonical")
    CP = mpatches.Patch(color="green", label="Canonical")
    BP = mpatches.Patch(color="black", label="Backbone")
    SP = mpatches.Patch(color="red", label="Stacking")
    # plt.legend(handles=[NCP,CP,BP,SP],prop={'size': 16})
    if savefig is None:
        plt.show()
    else:
        plt.tight_layout()
        plt.savefig(savefig)

    plt.clf()


if __name__ == '__main__':
    with open(PATH_CLUSTERS, 'rb') as f:
        data = pickle.load(f)

    try:
        os.mkdir(PATH_FIGURES)
    except FileExistsError:
        pass

    for i, c in enumerate(data):
        print(f"cluster {i} was found {len(c['names'])} times through all RNA 3D structures")
        path_to_draw = os.path.join(PATH_FIGURES, f'cluster_{i}.pdf')
        draw_graph(c['graph'], savefig=path_to_draw)

