import networkx as nx
import matplotlib.pyplot as plt
import graphviz
from graphviz import Digraph
# import numpy as np
# from mayavi import mlab
from networkx.drawing.nx_pydot import write_dot

def weightcheck(st):
    if st == 'Undirected Link':
        return 0.5
    elif st == 'Positive Increase':
        return 1
    elif st == 'Positive Decrease':
        return 1
    elif st == 'Directed Link':
        return 1
    else:
        pass
    return 0
def addpnode(entity):
    cnt = entity.count('&') + 1
    pts = entity.split('&')
    ptc = 'PX_'
    type = ''
    modi = ''
    iso = ''
    for g_cnt in range(cnt):
        g_items = pts[g_cnt].split('#')
        g_name = g_items[0]
        g_type = g_items[1]
        g_modi = g_items[2]
        g_iso = g_items[3]
        if g_cnt == 0:
            ptc = ptc + g_name
            type = type + g_type
            modi = modi + g_modi
            iso = iso + g_iso
        else:
            ptc = ptc + '&' + g_name
            type = type + '&' + g_type
            modi = modi + '&' + g_modi
            iso = iso + '&' + g_iso
    G.add_node(ptc)
    return ptc
def addgnode(entity):
    g_items = entity.split('#')
    g_name = g_items[0]
    g_type = g_items[1]
    g_modi = g_items[2]
    g_iso = g_items[3]
    # G.node(g_name, type=g_type, modification=g_modi, isoform=g_iso)
    G.add_node(g_name)
    return g_name
# main function
if __name__ == '__main__':
    G = nx.DiGraph()
    f = open("./input/all_level0.txt", encoding='utf-8')
    lines = f.readlines()
    beforeline = ''
    edgecnts = dict()
    gene_cnt = 0
    pt_cnt = 0
    pheno_cnt = 0
    bp_cnt = 0
    mf_cnt = 0
    cp_cnt = 0
    print('Start parsing.......')
    for line in lines:
        # 첫번째 줄 일때
        if 'Phenotype' in line:
            objs = line.split('\t')

            LHS_at = objs[0]
            LHS_entity = LHS_at[:-3]
            #protein complex일때
            if '&G' in LHS_entity:
                l_node = addpnode(LHS_entity)
            #gene일 때
            else:
                l_node = addgnode(LHS_entity)

            #association 확인
            asn = objs[1]
            edge_w = weightcheck(asn)


            RHS_at = objs[2]
            RHS_entity = RHS_at[:-3]

            #protein complex일때
            if '&G' in RHS_entity:
                r_node = addpnode(RHS_entity)
            #gene 혹은 function일 때
            else:
                r_node = addgnode(RHS_entity)

            contx = objs[3]
            ctx_li = contx.split('&')
            ctx_lis = ctx_li[0].split(':')

            #표현형(위암) PH id
            ph_id = ctx_lis[1]


            #association in_source
            ass_in = objs[4]
            #edge 추가
            if edge_w == 0.5:
                G.add_edge(l_node, r_node, weight='0.5', color='red')
                G.add_edge(r_node, l_node, weight='0.5', color='red')
            elif edge_w == -1:
                G.add_edge(l_node, r_node, weight='-1.0', color='red')
            else:
                G.add_edge(l_node, r_node, weight='1.0', color='red')
            if (l_node, r_node) in edgecnts:
                edgecnts[(l_node, r_node)] = edgecnts[(l_node, r_node)] + 1
            else:
                edgecnts[(l_node, r_node)] = 1

        # 두 번째 줄이면서 마지막 줄 일때
        elif 'Phenotype' in beforeline and 'Curation' in line:
            pass
        # 두 번째 줄이고 마지막이 아닐 때
        elif 'Phenotype' in beforeline:
            pass
        # 세 번째 줄일 때
        else:
            pass
        beforeline = line

    print('Finish parsing!')
    print('nodes are ' + str(G.number_of_nodes()))
    print('edges are ' + str(G.number_of_edges()))
    print('Start removing')
    H = G.copy()
    # for node in G.nodes():
    #     if not nx.has_path(H, 'PH00074992', node):
    #         H.remove_node(node)
    for key in edgecnts.keys():
        if edgecnts[key] < 5 and (key[0], key[1]) in H.edges():
            H.remove_edge(key[0], key[1])
    for node in G.nodes():
        if node == 'PH00074992':
            continue
        if not nx.has_path(H, 'PH00074992', node):
            H.remove_node(node)
        else:
            if nx.shortest_path_length(G, 'PH00074992', node) > 2:
                H.remove_node(node)
    print('Finish removing')
    for node in H.nodes():
        if 'PH' in node:
            pheno_cnt += 1
        elif 'PX' in node:
            pt_cnt += 1
        elif 'GE' in node:
            gene_cnt += 1
        elif 'CP' in node:
            cp_cnt += 1
        elif 'MF' in node:
            mf_cnt += 1
        elif 'BP' in node:
            bp_cnt += 1
        else:
            pass

    print('nodes are ' + str(H.number_of_nodes()))
    print('edges are ' + str(H.number_of_edges()))
    print('Phenotypes are ' + str(pheno_cnt) + ' Protein complex are ' + str(pt_cnt) + ' Genes are ' + str(gene_cnt) + ' Compounds are ' + str(cp_cnt)
    + ' MFs are ' + str(mf_cnt) + ' BPs are ' + str(bp_cnt))
    print('sum is ' + str(pheno_cnt+pt_cnt+gene_cnt+cp_cnt+ mf_cnt+bp_cnt))
    largest_hub = 'PH00074992'
    # Create ego graph of main hub
    hub_ego = nx.ego_graph(H, largest_hub)
    # Draw graph
    pos = nx.spring_layout(hub_ego)
    nx.draw(hub_ego, pos, node_color='b', edge_color='#C7C7CA', node_size=4, width = 0.5, arrowsize=0.5, with_labels=False)
    # Draw ego as large and red
    nx.draw_networkx_nodes(hub_ego, pos, nodelist=[largest_hub], node_size=100, node_color='r')
    # save_graph_as_svg(G, 'alpha_algorithm_dot')
    # nx.draw(G)
    plt.savefig("./output/networkGraph.png")
    # G.render('network_output2')
    # pos = nx.nx_agraph.graphviz_layout(H)
    # nx.draw(H, pos=pos)
    # write_dot(H, 'file.dot')
    # #3D graph
    # # reorder nodes from 0,len(G)-1
    # G=nx.convert_node_labels_to_integers(H)
    # # 3d spring layout
    # pos=nx.spring_layout(G,dim=3)
    # # numpy array of x,y,z positions in sorted node order
    # xyz=np.array([pos[v] for v in sorted(G)])
    # # scalar colors
    # scalars=np.array(G.nodes())+5
    #
    # mlab.figure(1, bgcolor=(0, 0, 0))
    # mlab.clf()
    #
    # pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
    #                     scalars,
    #                     scale_factor=0.1,
    #                     scale_mode='none',
    #                     colormap='Blues',
    #                     resolution=20)
    #
    # pts.mlab_source.dataset.lines = np.array(G.edges())
    # tube = mlab.pipeline.tube(pts, tube_radius=0.01)
    # mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
    #
    # mlab.savefig('mayavi2_spring.png')
