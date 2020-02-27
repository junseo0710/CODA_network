import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import graphviz
from graphviz import Digraph
from math import exp
# import numpy as np
# from mayavi import mlab
from networkx.drawing.nx_pydot import write_dot
disease_genes = ['PH00077705', 'PH00074992', 'PH00042592', 'PH00075982']
def graphst(node, st, graph):
    sum = 0
    for n in graph.nodes():
        if nx.has_path(graph, source = node, target = n) and nx.shortest_path_length(graph, source = node, target = n)==st:
            sum += 1
    return sum
def allst(node, dlist, graph):
    sum = 0
    for d in dlist:
        if nx.has_path(graph, source = node, target = d):
            sum += len(list(nx.all_shortest_paths(graph, source = node, target = d)))
    return sum
def Rad(node, dm, graph):
    sum = 0
    for g in graph.nodes():
        if nx.has_path(graph, source = node, target = g):
            sum += (dm + 1 - nx.shortest_path_length(graph, source= node, target = g))
        else:
            pass
    return sum
def LR(node, graph):
    sum = 0
    for dg in graph.nodes():
        if nx.has_path(graph, source = node, target = dg):
            sum += nx.shortest_path_length(graph, source = node, target = dg)
    return sum
def findcenter(graph, dlist):
    min = 99999
    for node in graph.nodes():
        sum = 0
        for dis in dlist:
            if nx.has_path(graph, source=node, target=dis):
                sum += nx.shortest_path_length(graph, source=node, target=dis)
        if sum < min:
            retnodelist = list()
            min = sum
            retnodelist.append(node)
        elif sum == min:
            retnodelist.append(node)
        else:
            continue
    return retnodelist
def stpathcheck(graph, node, dlist):
    minlen = 100
    sum = 0
    minst = list()
    for dis in dlist:
        if nx.has_path(graph, source=node, target=dis):
            spl = nx.shortest_path_length(graph, source=node, target=dis)
            if spl < minlen:
                minlen = spl
            sum += spl
            minst.append(minst)
        else:
            continue
    return sum/len(minst)
def stpathcheck2(graph, node, dlist, retnodelist):
    minlen = 100
    sum = 0
    cnt = 0
    minst = list()
    for cent in retnodelist:
        if nx.has_path(graph, source=cent, target=node):
            cnt += 1
            spl = nx.shortest_path_length(graph, source=cent, target=node)
            print(spl)
            if spl < minlen:
                minlen = spl
            sum += spl
        else:
            continue
    return minlen
def pathcheck(graph, node, dlist):
    for dis in dlist:
        if nx.has_path(graph, source=node, target=dis):
            return True
        else:
            continue
    return False
def weightcheck(st):
    if st == 'Undirected Link':
        return 0.5
    elif st == 'Positive Increase':
        return 1
    elif st == 'Positive Decrease':
        return -1
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
    G.add_node(g_name)
    return g_name


def targetchoose(graph, node):
    D = graph.degree[node]
    degreesum = 0
    for connode in graph.neighbors(node):
        degreesum = degreesum + 1/graph.degree[connode]
    return (1/D) * (1/degreesum)
def drugparsing(graph):
    ff = open("./input/basic_parsing_result.tsv", mode='r', encoding='utf-8')
    first = ff.readline()
    drug_lines = ff.readlines()
    ff.close()
    total_cnt = len(drug_lines)
    drug_dict = dict()
    target_score = dict()
    target_length = dict()
    target_degree = dict()
    drug_cnt = 0
    temp_cnt = 0
    rf = open("./output/drug_result_v8_final.txt", mode='w')
    highlist = ['DB00734', 'DB00751', 'DB01200', 'DB01267', 'DB04946', 'DB04948', 'DB06216', 'DB06288', 'DB08815', 'DB13988', 'DB00216', 'DB00247', 'DB00543', 'DB06148', 'DB09068', 'DB09225', 'DB00477', 'DB00421']

    target_dict = dict()
    str_dict = dict()
    retnodel = findcenter(graph, disease_genes)
    print(retnodel)
    # diameter = nx.diameter(graph)
    for drug in drug_lines:
        # print('++++++++++++++++++')
        temp_cnt += 1
        if temp_cnt%100 == 0:
            print("---------------------------------------")
            print(str(temp_cnt) + '/' + str(total_cnt))
            print("---------------------------------------")
        st_length = 0
        items = drug.split('\t')
        db_id = items[0]
        db_accuracy = items[3]
        # if 'approved' in db_accuracy:
        #     acc = 1
        # else:
        #     continue
        #     acc = 0
        target_id = items[5]
        target_cnt = target_id.count('|') + 1
        targetstring = ""
        if target_cnt == 1:
            target_id = target_id.split('\n')[0]
            if target_id == 'None':
                target_cnt -= 1
                continue
            elif target_id in graph.nodes() and graph.degree[target_id] < 5:
                continue
            elif target_id in target_dict:
                drug_dict[db_id] = target_dict[target_id]
                str_dict[db_id] = target_id
            elif target_id in graph.nodes() and pathcheck(graph, target_id, disease_genes):

                drug_cnt += 1
                # st_length = exp(stpathcheck2(graph, target_id, disease_genes, retnodel))**2
                real_length = stpathcheck2(graph, target_id, disease_genes, retnodel)
                st_length = exp(real_length)**2
                target_dict[target_id] = round((1 / st_length) * targetchoose(graph, target_id), 10)
                drug_dict[db_id] = target_dict[target_id]
                # drug_dict[db_id] = round((1 / st_length), 10)
                str_dict[db_id] = target_id
                target_dict[target_id] = drug_dict[db_id]
                target_score[target_id] = targetchoose(graph, target_id)
                target_degree[target_id] = graph.degree[target_id]
                target_length[target_id] = st_length

            else:
                target_cnt -= 1
                continue
        else:
            targets_list = target_id.split('|')
            targets_list[target_cnt-1] = targets_list[target_cnt-1].split('\n')[0]
            tot_val = 0
            temppath_val = 0
            target_cnt = 0
            for target_id in targets_list:

                if target_id == 'None':
                    continue
                elif target_id in graph.nodes() and graph.degree[target_id] < 5:
                    continue
                elif target_id in target_dict:
                    target_cnt += 1
                    tot_val = tot_val + target_dict[target_id]
                    if targetstring == "":
                        targetstring = target_id
                    else:
                        targetstring = targetstring + "+" + target_id
                elif target_id in graph.nodes() and pathcheck(graph, target_id, disease_genes):
                    target_cnt += 1
                    real_length = stpathcheck2(graph, target_id, disease_genes, retnodel)
                    st_length = exp(real_length)**2
                    target_dict[target_id] = round((1 / st_length) *  targetchoose(graph, target_id), 10)
                    # target_dict[target_id] = round((1 / st_length), 10)
                    if targetstring == "":
                        targetstring = target_id
                    else:
                        targetstring = targetstring + "+" + target_id
                    tot_val = tot_val + target_dict[target_id]
                    target_score[target_id] = targetchoose(graph, target_id)
                    target_degree[target_id] = graph.degree[target_id]
                    target_length[target_id] = st_length
                else:
                    continue
            if target_cnt == 0:
                continue
            # print("total value is " + str(tot_val))
            drug_dict[db_id] = round((tot_val/target_cnt), 10)
            str_dict[db_id] = targetstring
        if target_cnt == 0:
            continue
        # print(drug_dict[db_id])
        rf.write(db_id + '\t' + str(drug_dict[db_id]) + '\t' + str_dict[db_id] + '\n')
    rf.close()
    tf = open("./output/target_result_v8_final.txt", mode='w')
    for key in target_score.keys():
        tf.write(key + '\t' + str(target_length[key]) + '\t' + str(target_score[key]) + '\t' + str(target_degree[key]) + '\n')
    tf.close()
# main function
if __name__ == '__main__':
    G = nx.Graph()
    f = open("./input/stomach_neoplasm_level1.txt", encoding='utf-8')
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
                G.add_edge(l_node, r_node, weight=1.0, color='red')
                G.add_edge(r_node, l_node, weight=1.0, color='red')
            elif edge_w == -1:
                G.add_edge(l_node, r_node, weight=1.0, color='red')
            else:
                G.add_edge(l_node, r_node, weight=1.0, color='red')
            # if set([l_node, r_node]) in edgecnts:
            #     edgecnts[set([l_node, r_node])] = edgecnts[set([l_node, r_node])] + 1
            # else:
            #     edgecnts[set([l_node, r_node])] = 1

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
    for node in G.nodes():
        if node in disease_genes:
            continue
        if not node in H.nodes():
            continue
        if 'PH' in node:
            H.remove_node(node)
            continue
    # for key in edgecnts.keys():
    #     if edgecnts[key] < 2 and (key[0], key[1]) in H.edges():
    #         H.remove_edge(key[0], key[1])
    R = H.copy()
    for node in R.nodes():
        if node in disease_genes:
            continue
        if not node in H.nodes():
            continue
        if nx.has_path(H, node, disease_genes[0]) or nx.has_path(H, node, disease_genes[1]) or nx.has_path(H, node, disease_genes[2]) or nx.has_path(H, node, disease_genes[3]):
            if stpathcheck(H, node, disease_genes) > 3:
                H.remove_node(node)
            else:

                continue
            # if nx.shortest_path_length(H, source=node, target=disease_gene) > 4:
            #     H.remove_node(node)
            continue
        else:
            H.remove_node(node)
    # K = H.copy()
    # for node in K.nodes():
    #     if node in disease_genes:
    #         continue
    #     if not node in H.nodes():
    #         continue
    #     if not 'GE' in node:
    #         H.remove_node(node)
    #     else:
    #         continue
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
    pos3 = nx.spring_layout(G)
    nx.draw(G, pos3, node_color='b', edge_color='#C7C7CA', node_size=4, width = 0.5, arrowsize=0.5, with_labels=False)
    nx.draw_networkx_nodes(G, pos3, nodelist=['GE05350394', 'GE05347819', 'GE05350268', 'GE05351705', 'GE05350265', 'GE05361827', 'GE05347980', 'GE05344450', 'GE05350143', 'GE05345033', 'GE05346235'],
     node_size=10, node_color='g')
    # nx.draw_networkx_nodes(H, pos3, nodelist=centernodes, node_size=10, node_color='#BFAFF1')
    nx.draw_networkx_nodes(G, pos3, nodelist=[disease_genes[0]], node_size=100, node_color='r', node_shape='^')
    nx.draw_networkx_nodes(G, pos3, nodelist=[disease_genes[1]], node_size=100, node_color='r')
    nx.draw_networkx_nodes(G, pos3, nodelist=[disease_genes[2]], node_size=100, node_color='r', node_shape='s')
    nx.draw_networkx_nodes(G, pos3, nodelist=[disease_genes[3]], node_size=100, node_color='r', node_shape='d')
    plt.savefig("./output/networkGraph8_fake.png")
    # drugparsing(H)
