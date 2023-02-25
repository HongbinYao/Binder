from copy import deepcopy

class Vertex:
    def __init__(self, ver):
        self.ver_inf = ver
        self.value = 0
        self.adjacency = []
        self.edge_value = {}
    
    def addneighbor(self, nbr):#添加邻接点
        if nbr not in self.edge_value:
            self.edge_value[nbr] = 1
        else:
            self.edge_value[nbr] += 1

    def __str__(self):
        return str(self.ver_inf) + 'adjacencyTo' + str([x.ver_inf for x in self.adjacency])

    def get_adjacency(self):#返回该点的所有邻接点
        return self.adjacency
    
    def get_weight(self, nbr):#返回到该邻接点的边的权重
        return self.edge_value[nbr]
    
    def get_edge(self):
        for adjacency in self.edge_value:
            print((self.ver_inf,adjacency),'-->',self.edge_value[adjacency],end=' ')
        print('')
        return None

    
class adjacency_list:
    def __init__(self):
        self.vers_class = {}
        self.vers_inf = []
        self.vers_num = 0
        self.end_point = {}

    def addVertex(self,ver):#只修改点的权重
        if ver not in self.vers_inf:
            self.vers_num = self.vers_num + 1
            self.vers_inf.append(ver)
            self.vers_class[ver] = Vertex(ver); self.vers_class[ver].value += 1
        else:
            self.vers_class[ver].value += 1

    def addEdge(self,front,back):#只修改边的权重
        if front not in self.vers_inf:
            self.addVertex(front)
            self.vers_class[front] = Vertex(front)
            self.vers_class[front].addneighbor(back)#修改边权重
            self.vers_class[front].adjacency.append(back)
        else:
            self.vers_class[front].addneighbor(back)#修改边权重
            if back not in self.vers_class[front].adjacency:
                self.vers_class[front].adjacency.append(back)
        if back not in self.vers_inf:
            self.vers_inf.append(back)

    def get_end_point(self,left,right):
        temp = []
        if left not in self.end_point:
            temp.append(right)
            self.end_point[left] = temp
        else:
            if right not in self.end_point[left]:
                temp.append(right)
                self.end_point[left] += temp

    def get_vers_inf(self):
        return self.vers_inf
    
    
def construct_graph(assembly_exon): #将genome_guided_exon添加到splice graph中
    splice_graph = adjacency_list()

    for tran in assembly_exon:
        for index in range(0,len(tran)-1):
            front = tran[index]; back = tran[index+1]
            splice_graph.addVertex(front); splice_graph.addVertex(back)
            splice_graph.addEdge(front, back)
    
        splice_graph.get_end_point(tran[0],tran[-1])
    
    return splice_graph


def rectify(de_nove_temp, splice_graph):#rectify de_nove exon
    #提取splice graph中的genome_guided exon 点
    genome_guided_exon = deepcopy(splice_graph.vers_inf)
    genome_guided_exon.sort()
    
    #提取genome_guided 剪切位点
    genome_guided_junction = []
    for exon in genome_guided_exon:
        genome_guided_junction.append(exon[0]); genome_guided_junction.append(exon[1])
    genome_guided_junction = list(set(genome_guided_junction)), genome_guided_junction.sort()
    
    #转换de_nove exon 格式
    de_nove_exon = []
    for tran in de_nove_temp:
        temp = []
        for index in range(0,len(tran)-2, 2):
            exon = ()
            exon += (tran[index],); exon += (tran[index + 1],)
            temp.append(exon)
        de_nove_exon.append(temp)
        
    #对de_nove exon 进行修正
    de_nove_rectify = []
    not_append_end = []
    for tran in de_nove_exon:
        if tran[0][1] not in genome_guided_junction and tran[-1][0] not in genome_guided_junction:
            not_append_end.append(de_nove_exon.index(tran))
        temp = []
        for exon in tran:
            if exon in genome_guided_exon:#若exon在genome_guided 中,直接添加该exon
                temp.append(exon)
                continue
            exon_index = 0
            while exon[0] >= genome_guided_exon[exon_index][1]:
                exon_index += 1
                #终止条件：已经遍历到genome_guided 最后一个exon
                if exon_index >= len(genome_guided_exon):
                    #此时说明de_nove exon 遍历到genome_guided 最后一个exon 依然在右侧 则不添加该双端exon
                    not_append_end.append(de_nove_exon.index(tran))
                    break
            else:
                while exon[1] > genome_guided_exon[exon_index][1]:
                    exon_index += 1
                    #终止条件：已经遍历到genome_guided 最后一个exon
                    if exon_index >= len(genome_guided_exon):
                        if genome_guided_exon[exon_index-1][0] < tran[tran.index(exon)-1][1]:
                            #此时de_nove最后一个exon要修正为的exon 与de_nove前一个exon相交
                            # temp.append(exon)
                            break
                        else:
                            #此时de_nove exon的左边界 >= genome_guided exon的左边界 直接添加即可
                            if exon[0] >= genome_guided_exon[exon_index-1][0]:
                                temp.append(genome_guided_exon[exon_index-1])
                                break
                            #此时de_nove exon的左边界 < genome_guided exon的左边界 应继续向前探索
                            last_exon_index = exon_index
                            while exon[0] <= genome_guided_exon[exon_index-1][0]:
                                exon_index -= 1
                            else:
                                for index in range(exon_index, last_exon_index):
                                    temp.append(genome_guided_exon[index])
                                break
                else:#exon[1] <= genome_guided_exon[exon_index][1]
                    if exon[1] > genome_guided_exon[exon_index][0]:
                        # if exon[0] < genome_guided_exon[exon_index][0]:
                        #     #此时de_nove exon的左边界 < genome_guided exon的左边界 应继续向前探索
                        #     if exon_index == 0:#此时直接添加第一个genome_guided exon即可
                        #         temp.append(genome_guided_exon[exon_index])
                        #     #此时de_nove exon左边界与前一个genome_guided exon左边界相同
                        #     #de_nove exon右边界与当前genome_guided exon右边界相同 且de_nove exon值较小时
                        #     elif exon[0] == genome_guided_exon[exon_index-1][0] and exon[1] == genome_guided_exon[exon_index][1] and (exon[1] - exon[0]) < 301:#634
                        #         temp.append(exon)
                        #     #此时de_nove exon左边界 <= genome_guided exon左边界 可以继续向前探索
                        #     elif exon[0] <= genome_guided_exon[exon_index-1][0]:
                        #         last_exon_index = exon_index
                        #         while exon[0] <= genome_guided_exon[exon_index][0] and exon_index >= 0:
                        #             # if exon_index == 0:
                        #             #     break
                        #             exon_index -= 1
                        #         else:
                        #             # if exon_index == 0:
                        #             #     temp.append(genome_guided_exon[exon_index])
                        #             exon_index += 1
                        #             while exon_index <= last_exon_index:
                        #                 #若前后两个genome_guided exon有交集
                        #                 if exon_index+1 <= last_exon_index and genome_guided_exon[exon_index][1] >= genome_guided_exon[exon_index+1][0]:
                        #                     #判断de_nove exon 与 genome_guided exon
                        #                     if exon[0] == genome_guided_exon[exon_index][0] or exon[1] == genome_guided_exon[exon_index][1]:
                        #                         temp.append(genome_guided_exon[exon_index]); exon_index += 2
                        #                     else:#exon[1] == genome_guided_exon[exon_index+1][1]
                        #                         temp.append(genome_guided_exon[exon_index+1]); exon_index += 2
                        #                 else:
                        #                     temp.append(genome_guided_exon[exon_index]); exon_index += 1
                        #     else:#exon[0] >= genome_guided_exon[exon_index-1][1]
                        #         temp.append(genome_guided_exon[exon_index])
                        # elif (genome_guided_exon[exon_index][1]-genome_guided_exon[exon_index][0]) > 2*(exon[1]-exon[0]) and 1 <= tran.index(exon) <= len(tran) - 3 and (exon[1]-exon[0]) > 54:
                        #     #如果要修正为的genome_guided 的exon与de_nove exon相差过大 那么添加de_nove exon
                        #     temp.append(exon)
                        # elif genome_guided_exon[exon_index][0] < tran[tran.index(exon)-1][1] and 1 <= tran.index(exon) <= len(tran) - 2:
                        #     #如果要修正为genome_guided 的exon与前一个de_nove exon相交 此时添加de_nove exon
                        #     temp.append(exon)
                        # elif genome_guided_exon[exon_index][0] < tran[tran.index(exon)-1][1] and tran.index(exon) == len(tran) - 1:
                        #     #如果要修正为genome_guided 的exon与前一个de_nove exon相交 且为最后一个de_nove exon 那么pass
                        #     continue
                        # elif 1 <= tran.index(exon) <= len(tran) - 2 and genome_guided_exon[exon_index][1] >= tran[tran.index(exon)+1][0]:
                        #     #若要修正为genome_guided exon与de_nove 后一个exon相交 那么保留de_nove exon或删除
                        #     if exon_index == len(genome_guided_exon) - 1:#此时已经遍历到最后一个genome_guided exon
                        #         # if exon[0] == genome_guided_exon[exon_index][0]:
                        #         temp.append(genome_guided_exon[exon_index])
                        #     else:
                        #         temp.append(exon)
                        # else:#其他情况
                        #     #例如 exon[0] >= genome_guided_exon[exon_index][0]
                        #     temp.append(genome_guided_exon[exon_index])
                            
                        
                        while exon[1] == genome_guided_exon[exon_index][1]:#右边界相同
                            if exon[0] == genome_guided_exon[exon_index-1][0]:
                                temp.append(exon)
                                break
                            else:
                                exon_index += 1
                                if exon_index >= len(genome_guided_exon):
                                    temp.append(genome_guided_exon[exon_index-1])
                                    break
                        else:
                            if exon[1] <= genome_guided_exon[exon_index][0]:#迭代到最后没有与de_nove exon相同的exon
                                if genome_guided_exon[exon_index-1][0] < tran[tran.index(exon)-1][1] and 0 < tran.index(exon) <= len(tran) - 2:
                                    #如果要修正为genome_guided 的exon与前一个de_nove exon相交 此时添加de_nove exon
                                    temp.append(exon)
                                elif exon[0] <= genome_guided_exon[exon_index-2][0] and 0 < tran.index(exon) <= len(tran) - 2 and exon_index >= 3:
                                    #此时exon不光与当前genome_guided exon相交还与前一个genome_guided exon相交
                                    #而且当前genome_guided exon不与前一个de_nove exon相交
                                    #因为当前de_nove exon与前一个genome_guided exon左边界相交 说明前一个de_nove exon与前一个genome_guided exon不相交
                                    if genome_guided_exon[exon_index-1][0] <= genome_guided_exon[exon_index-2][1]:
                                        #此时当前genome_guided exon与前一个genome_guided exon 相交 应当添加前前个genome_guided exon
                                        temp.append(genome_guided_exon[exon_index-3])
                                        temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                                        continue
                                    else:#此时当前genome_guided exon与前一个genome_guided exon 不相交
                                        temp.append(genome_guided_exon[exon_index-2])
                                        temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                                        continue
                                else:
                                    temp.append(genome_guided_exon[exon_index-1])
                                    continue
                            elif genome_guided_exon[exon_index][0] < tran[tran.index(exon)-1][1] and 0 < tran.index(exon) <= len(tran) - 2:
                                #如果要修正为genome_guided 的exon与前一个de_nove exon相交 此时添加de_nove exon
                                temp.append(exon)
                            elif genome_guided_exon[exon_index][0] < tran[tran.index(exon)-1][1] and tran.index(exon) == len(tran) - 1:
                                #如果要修正为genome_guided 的exon与前一个de_nove exon相交 且为最后一个de_nove exon 那么pass
                                continue
                            elif (genome_guided_exon[exon_index][1]-genome_guided_exon[exon_index][0]) > 2*(exon[1]-exon[0]) and 0 < tran.index(exon) <= len(tran) - 3 and (exon[1]-exon[0]) > 54:
                                #如果要修正为的genome_guided 的exon与de_nove exon相差过大 那么添加de_nove exon
                                temp.append(exon)
                            elif exon[0] <= genome_guided_exon[exon_index-1][0] and 0 < tran.index(exon) <= len(tran) - 2 and exon_index >= 2:
                                #此时exon不光与当前genome_guided exon相交还与前一个genome_guided exon相交
                                #而且当前genome_guided exon不与前一个de_nove exon相交
                                #因为当前de_nove exon与前一个genome_guided exon左边界相交 说明前一个de_nove exon与前一个genome_guided exon不相交
                                if genome_guided_exon[exon_index][0] < genome_guided_exon[exon_index-1][1]:
                                    #此时当前genome_guided exon与前一个genome_guided exon 相交 应当添加前前个genome_guided exon
                                    temp.append(genome_guided_exon[exon_index-2])
                                    temp.append((genome_guided_exon[exon_index][0],exon[1]))
                                    continue
                                else:#此时当前genome_guided exon与前一个genome_guided exon 不相交
                                    temp.append(genome_guided_exon[exon_index-1])
                                    temp.append((genome_guided_exon[exon_index][0],exon[1]))
                                    continue
                            elif 0 < tran.index(exon) <= len(tran) - 2 and genome_guided_exon[exon_index][1] >= tran[tran.index(exon)+1][0]:
                                #若要修正为genome_guided exon与de_nove 后一个exon相交 那么保留de_nove exon或删除
                                if exon_index >= len(genome_guided_exon) - 1:#此时已经遍历到最后一个genome_guided exon
                                    if exon[0] == genome_guided_exon[exon_index][0]:
                                        temp.append(genome_guided_exon[exon_index])
                                    else:
                                        temp.append(exon)
                                else:
                                    temp.append(exon)
                            else:#遍历到最后一个exon
                                # if tran.index(exon) > 0:
                                temp.append(genome_guided_exon[exon_index])
                    else:#exon[1] <= genome_guided_exon[exon_index][0]:
                        # if exon_index == 0:
                        # #此时说明exon的右边界 <= genome_guided 第一个exon的左边界
                        #     not_append_end.append(de_nove_exon.index(tran))
                        # elif exon[0] < genome_guided_exon[exon_index-1][0]:
                        #     if exon[0] == genome_guided_exon[exon_index-2][0] and exon[1] == genome_guided_exon[exon_index-1][1] and (exon[1] - exon[0]) < 301:#634
                        #         temp.append(exon)
                        #     elif exon[0] <= genome_guided_exon[exon_index-2][0]:
                        #         last_exon_index = exon_index-1
                        #         while exon[0] <= genome_guided_exon[exon_index-1][0] and exon_index-1 >= 0:
                        #             exon_index -= 1
                        #         else:
                        #             exon_index += 1
                        #             while exon_index <= last_exon_index:
                        #                 #若前后两个genome_guided exon有交集
                        #                 if exon_index+1 <= last_exon_index and genome_guided_exon[exon_index][1] >= genome_guided_exon[exon_index+1][0]:
                        #                     #判断de_nove exon 与 genome_guided exon
                        #                     if exon[0] == genome_guided_exon[exon_index][0] or exon[1] == genome_guided_exon[exon_index][1]:
                        #                         temp.append(genome_guided_exon[exon_index]); exon_index += 2
                        #                     else:#exon[1] == genome_guided_exon[exon_index+1][1]
                        #                         temp.append(genome_guided_exon[exon_index+1]); exon_index += 2
                        #                 else:
                        #                     temp.append(genome_guided_exon[exon_index]); exon_index += 1
                        #     else:#genome_guided_exon[exon_index-2][1] < exon[0] < genome_guided_exon[exon_index-1][0]
                        #         temp.append(genome_guided_exon[exon_index-1])
                        # elif genome_guided_exon[exon_index-2][1] <= exon[0] < genome_guided_exon[exon_index-1][0]:
                        #     if tran[tran.index(exon)-1][0] < genome_guided_exon[exon_index-1][0] < tran[tran.index(exon)-1][1]:
                        #         #如果要修正为genome_guided 的exon与前一个de_nove exon相交 此时添加de_nove exon
                        #         temp.append(exon)
                        #     elif genome_guided_exon[exon_index-1][0] < tran[tran.index(exon)-1][1] and tran.index(exon) == len(tran) - 1:
                        #         #如果要修正为genome_guided 的exon与前一个de_nove exon相交 且为最后一个de_nove exon 那么pass
                        #         continue
                        #     elif (genome_guided_exon[exon_index-1][1]-genome_guided_exon[exon_index-1][0]) > 2*(exon[1]-exon[0]) and 1 <= tran.index(exon) <= len(tran) - 3 and (exon[1]-exon[0]) > 54:
                        #         #如果要修正为的genome_guided 的exon与de_nove exon相差过大 那么添加de_nove exon
                        #         temp.append(exon)
                        #     else:
                        #         temp.append(genome_guided_exon[exon_index-1])
                        #其他情况 pass
                        
                        if exon_index == 0:#tran.index(exon) == 0
                            not_append_end.append(de_nove_exon.index(tran))
                        elif exon[0] <= genome_guided_exon[exon_index-1][1]:
                            if genome_guided_exon[exon_index-1][0] < tran[tran.index(exon)-1][1]:
                                #此时de_nove要修正为genome_guided exon与de_nove前一个exon相交
                                #not_append_end.append(de_nove_exon.index(tran))
                                temp.append(exon)
                            elif exon[0] <= genome_guided_exon[exon_index-2][0]:#如果de_nove exon横跨两个genome_guided exon 那么保留这两个exon
                                    if genome_guided_exon[exon_index-1][0] < genome_guided_exon[exon_index-2][1]:
                                        #此时当前genome_guided exon与前一个genome_guided exon 相交 应当添加前前个genome_guided exon
                                        if exon[0] <= genome_guided_exon[exon_index-3][0]:
                                            temp.append(genome_guided_exon[exon_index-3])
                                            temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                                        else:
                                            temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                                    else:
                                        temp.append(genome_guided_exon[exon_index-2])
                                        temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                            elif exon[0] > genome_guided_exon[exon_index-2][1]:#de_nove exon 只横跨一个exon
                                    if tran.index(exon) == len(tran) - 1:
                                        temp.append(genome_guided_exon[exon_index-1])
                                    else:
                                        if (exon[1]-exon[0]) - (genome_guided_exon[exon_index-1][1]-genome_guided_exon[exon_index-1][0]) > 0.25*(genome_guided_exon[exon_index-1][1]-genome_guided_exon[exon_index-1][0]):
                                            #如果de_nove exon的右边界超过1/4 genome_guided exon 那么添加de_nove exon的右边界
                                            temp.append((genome_guided_exon[exon_index-1][0],exon[1]))
                                        else:
                                            temp.append(genome_guided_exon[exon_index-1])
                            else:
                                temp.append(genome_guided_exon[exon_index-1])
                        else:
                            not_append_end.append(de_nove_exon.index(tran))
        temp = list(set(temp)); temp.sort()
        if len(temp) <= 2:
            continue
        de_nove_rectify.append(temp)
        
    not_append_end = list(set(not_append_end))
    return de_nove_rectify, not_append_end

        
def append_exon(de_nove_rectify, not_append_end, splice_graph):
    genome_guided_exon = deepcopy(splice_graph.vers_inf)
    genome_guided_exon.sort()
    
    pair_end = splice_graph.end_point
    temp = []
    for start, end in pair_end.items():
        for exon in end:
            temp.append(exon)
    
    for tran in de_nove_rectify:
        # if de_nove_rectify.index(tran) not in not_append_end:
        #     left = tran[0]; right = tran[-1]
        #     # if left in genome_guided_exon and right in genome_guided_exon and right in temp and (right[1] - left[0]) < 10*left[0] and (genome_guided_exon[-1][1] - genome_guided_exon[0][0]) > 4355:
        #         #当3' 5'端距离过大 #或3' 5'端不在genome_guided exon #或5'端不在pair_end中 #或genome_guided3' 5'端距离过小时不添加de_nove end_point
        #     splice_graph.get_end_point(left,right)
        #     for index in range(len(tran)-1):  
        #         front = tran[index]; back = tran[index + 1]
        #         splice_graph.addVertex(front); splice_graph.addVertex(back)
        #         # if back in genome_guided_exon and front in genome_guided_exon:
        #         #     if genome_guided_exon.index(back) - genome_guided_exon.index(front) <= 3:
        #         #         splice_graph.addEdge(front, back)
        #         # else:
        #         #     splice_graph.addEdge(front, back)
        #         splice_graph.addEdge(front, back)
        # else:
        #     for index in range(len(tran)-1):
        #         front = tran[index]; back = tran[index + 1]
        #         splice_graph.addVertex(front); splice_graph.addVertex(back)
        #         # if back in genome_guided_exon and front in genome_guided_exon:
        #         #     if genome_guided_exon.index(back) - genome_guided_exon.index(front) <= 3:
        #         #         splice_graph.addEdge(front, back)
        #         # else:
        #         #     splice_graph.addEdge(front, back)
        #         splice_graph.addEdge(front, back)
        if de_nove_rectify.index(tran) not in not_append_end:
            left = tran[0]; right = tran[-1]
            splice_graph.get_end_point(left,right)
            for index in range(len(tran)-1):
                front = tran[index]; back = tran[index + 1]
                splice_graph.addVertex(front); splice_graph.addVertex(back)
                splice_graph.addEdge(front, back)
        else:
            for index in range(len(tran)-1):
                front = tran[index]; back = tran[index + 1]
                splice_graph.addVertex(front); splice_graph.addVertex(back)
                splice_graph.addEdge(front, back)
        
    return splice_graph

if __name__ == '__main__':
    g = construct_graph(genome_guided_temp)#将genome_guided_exon添加到splice graph中
    de_nove_rectify = rectify(de_nove_temp, g)#对de_nove exon 进行修正
    g = append_de_nove_exon(de_nove_rectify, g)#添加修正后的de_nove exon
    path = path_search(g)
# graph = {}
# for i in g.vers_inf:
#     graph[i] = g.vers_class[g.vers_inf.index(i)].adjacency                 
# for i in g.vers_inf:
#     print(i,g.vers_class[g.vers_inf.index(i)].adjacency,'\n')
# print(g.vers_inf)
# print(g.end_point)
#计算可变剪接的种类
