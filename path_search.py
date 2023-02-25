# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:05:39 2022

@author: yhb
"""
from copy import deepcopy
from adjacency_list import *

def dfs(vertice,end,graph,trace):
    trace = deepcopy(trace)# 深拷贝，对不同起点，走过的路径不同       
    
    if vertice not in trace:
        trace.append(vertice)
        if vertice == end:
            path.append(trace)
            return None
    
    # 深度优先递归遍历
    for successor_node in graph[vertice]:
        dfs(successor_node,end,graph,trace)


def recursion(splice_graph):  
    path = []
    ends = splice_graph.end_point
    
    graph = {}
    for i in splice_graph.vers_inf:
        graph[i] = splice_graph.vers_class[splice_graph.vers_inf.index(i)].adjacency
    
    for k in splice_graph.vers_inf:
        if k in ends:
            for end in ends[k]:
                dfs(k,end,graph,[])
    print(path)
    return path


class Stack:
    def __init__(self):
        self.stack = []
        
    def push(self, element):
        self.stack.append(element)
        
    def pop(self):
        return self.stack.pop()
            
    def get_top(self):
        if len(self.stack) > 0:
            return self.stack[-1]
        else:
            return None
        
# def simply_path_search(splice_graph):
#     alternative_exon = []
#     alternative_edge = []
#     alternative_front = []
#     # splice_num = 0
#     graph = {}
#     for i in splice_graph.vers_inf:
#         graph[i] = splice_graph.vers_class[splice_graph.vers_inf.index(i)].adjacency
        
#     graph = {k:sorted(v) for k, v in graph.items()}
#     graph = dict(sorted(graph.items(),key=lambda x:x[0]))
#     # front_exon = {}
#     for front, back in graph.items():
#         back.sort()
#         if len(back) >= 2:
#             alternative_front.append(front)
#             # splice_num += 1
#             if len(back) == 2:
#                 if max(back[0][0],back[1][0]) < min(back[0][1],back[1][1]):
#                     alternative_exon.append(back[0]); alternative_exon.append(back[-1])
#                 else:
#                     alternative_edge.append(back[-1]); alternative_exon.append(back[0])
#                     # front_exon[back[-1]] = front
#             else:#len(back) >= 3
#                 deepcopy = deepcopy(back); deepcopy.sort()
#                 temp = deepcopy.pop(0)#把第一个可变剪接赋值给temp 由于对deepcopy经过了排序 所以temp必定是可变剪接点
#                 while deepcopy:
#                     if max(deepcopy[0][0],temp[0]) < min(deepcopy[0][1],temp[1]):
#                         #如果此时temp和下一个可变剪接相交 说明两个都为可变剪接点
#                         alternative_exon.append(temp); alternative_exon.append(deepcopy[0])
#                         temp = deepcopy.pop(0)#此时把temp赋值为下一个可变剪接 用于对下下个可变剪接做判断
#                     else:#此时temp为可变剪接点 下一个为可变剪接边 temp保持不变
#                         alternative_edge.append(deepcopy[0]); alternative_exon.append(temp)
#                         # front_exon[deepcopy[0]] = front
#                         deepcopy.pop(0)
    
#     chief_stack = deepcopy(Stack())
#     deputy_stack = deepcopy(Stack())
#     graph_exon = list(graph.keys()); graph_exon.sort()
#     ends = splice_graph.end_point
#     path = []
    
#     for vertice in graph_exon:
#         if vertice in ends:
#             for end in ends[vertice]:
#                 anthor = False
#                 alternative_distance = []#用于计算两个可变剪接边之间的距离
#                 alternative_distance.append(vertice)
#                 while alternative_distance[-1] != end:
#                     if alternative_distance[-1] in alternative_front:
#                         back = graph[alternative_distance[-1]]
#                         if len(back) == 2:
#                             if max(back[0][0],back[1][0]) < min(back[0][1],back[1][1]):
#                                 if anthor == False:
#                                     alternative_distance.append(back[0])#; alternative_distance.append(back[1])
#                                 else:
#                                     alternative_distance.append(back[1]); anthor = False
#                             else:#此时只添加可变剪接边
#                                 if anthor == False:
#                                     alternative_distance.append(back[-1])
#                                 else:
#                                     alternative_distance.append(back[0]); anthor = False
#                         else:
#                             deepcopy = deepcopy(back); deepcopy.sort()
#                             temp = deepcopy.pop(0)#把第一个可变剪接赋值给temp 由于对deepcopy经过了排序 所以temp必定是可变剪接点
#                             while deepcopy:
#                                 if max(deepcopy[0][0],temp[0]) < min(deepcopy[0][1],temp[1]):
#                                     #如果此时temp和下一个可变剪接相交 说明两个都为可变剪接点
#                                     if anthor == False:
#                                         alternative_distance.append(deepcopy(temp))#; alternative_exon.append(deepcopy[0])
#                                         temp = deepcopy.pop(0)#此时把temp赋值为下一个可变剪接 用于对下下个可变剪接做判断
#                                     else:
#                                         alternative_distance.append(deepcopy[0])
#                                         temp = deepcopy.pop(0); anthor = False
#                                 else:#此时temp为可变剪接点 下一个为可变剪接边 temp保持不变 此时只添加可变剪接边
#                                     if anthor == False:
#                                         alternative_distance.append(deepcopy.pop(0))#; alternative_exon.append(temp)
#                                     else:
#                                         alternative_distance.append(temp); anthor = False
#                     else:
#                         if graph[alternative_distance[-1]]:
#                             alternative_distance.append(graph[alternative_distance[-1]][0])
#                         else:#说明没有到达end 但已经没有路了 此时需要回退
#                             alternative_distance.pop()
#                             while graph[alternative_distance[-1]][0] > end:
#                                 alternative_distance.pop()
#                             else:#此时说明当前最后一个点小于end
#                                 anthor = True
                         
#                 append = False
#                 temp = []
#                 chief_stack.push(vertice); deputy_stack.push(graph[vertice])
#                 while chief_stack.stack:
#                     adjacency_vertices = deepcopy(deputy_stack.pop())
#                     if adjacency_vertices:
#                         if adjacency_vertices[0] in alternative_edge and chief_stack.get_top() in alternative_front and adjacency_vertices[0] <= end:
#                             #判断要添加的exon 是否和前一个 exon 组成可变剪接边
#                             # last_front = chief_stack.get_top()
#                             if append == False:#此时前面没有可变剪接边
#                                 temp.append(adjacency_vertices[0])
#                                 append = True
#                             else:#append == True:
#                                 if chief_stack.get_top() < temp[-1]:
#                                     continue
#                                 else:
#                                     #此时前面已经有一个可变剪接边 需要判断此时是否能添加另一个可变剪接边
#                                     if alternative_distance.index(adjacency_vertices[0]) - alternative_distance.index(temp[-1]) >= 5:
#                                         #此时两个可变剪接边的距离大于等于10 则后一个可变剪接边不可以添加
#                                         adjacency_vertices.pop(0)
#                                     elif alternative_distance.index(adjacency_vertices[0]) - alternative_distance.index(temp[-1]) <= 4:
#                                         #此时两个可变剪接边的距离小于等于3 则可以同时存在
#                                         temp.append(adjacency_vertices[0])
#                         elif adjacency_vertices[0] > end:
#                             adjacency_vertices.pop(0)
                            
#                         if adjacency_vertices:
#                             chief_stack.push(adjacency_vertices.pop(0)); deputy_stack.push(adjacency_vertices)
#                             deputy_stack.push(graph[chief_stack.get_top()])
#                         else:
#                             deputy_stack.push(adjacency_vertices)
#                             if chief_stack.get_top() == end:
#                                 # print(chief_stack.get_top());print(vertice)
#                                 path.append(deepcopy(chief_stack.stack))
#                                 chief_stack.pop(); deputy_stack.pop()
#                                 # append = False
#                             else:
#                                 if chief_stack.get_top() in temp:
#                                     temp.remove(chief_stack.get_top())
#                                 if not temp:
#                                     append = False
#                                 chief_stack.pop(); deputy_stack.pop()
                        
#                     else:
#                         deputy_stack.push(adjacency_vertices)
#                         if chief_stack.get_top() == end:
#                             # print(chief_stack.get_top());print(vertice)
#                             path.append(deepcopy(chief_stack.stack))
#                             chief_stack.pop(); deputy_stack.pop()
#                             # append = False
#                         else:
#                             if chief_stack.get_top() in temp:
#                                 temp.remove(chief_stack.get_top())
#                             if not temp:
#                                 append = False
#                             chief_stack.pop(); deputy_stack.pop()
    
#     return path
def simply_path_search(splice_graph):
    alternative_exon = []
    alternative_edge = []
    alternative_front = []
    
    graph = {}
    for i in splice_graph.vers_inf:
        graph[i] = splice_graph.vers_class[splice_graph.vers_inf.index(i)].adjacency
        
    graph = {k:sorted(v,reverse=True) for k, v in graph.items()}; graph = dict(sorted(graph.items(),key=lambda x:x[0]))
    
    for front, back in graph.items():
        back.sort()
        if len(back) >= 2:
            alternative_front.append(front)
            if len(back) == 2:
                if max(back[0][0],back[1][0]) < min(back[0][1],back[1][1]):
                    alternative_exon.append(back[0]); alternative_exon.append(back[-1])
                else:
                    alternative_edge.append(back[-1]); alternative_exon.append(back[0])
            else:#len(back) >= 3
                copy = deepcopy(back); copy.sort()
                temp = copy.pop(0)#把第一个可变剪接赋值给temp 由于对deepcopy经过了排序 所以temp必定是可变剪接点
                while copy:
                    if max(copy[0][0],temp[0]) < min(copy[0][1],temp[1]):
                        #如果此时temp和下一个可变剪接相交 说明两个都为可变剪接点
                        alternative_exon.append(temp); alternative_exon.append(copy[0])
                        temp = copy.pop(0)#此时把temp赋值为下一个可变剪接 用于对下下个可变剪接做判断
                    else:#此时temp为可变剪接点 下一个为可变剪接边 temp保持不变
                        alternative_edge.append(copy[0]); alternative_exon.append(temp)
                        copy.pop(0)
    
    chief_stack = deepcopy(Stack())
    deputy_stack = deepcopy(Stack())
    graph_exon = list(graph.keys()); graph_exon.sort()
    ends = splice_graph.end_point; ends = dict(sorted(ends.items(),key=lambda x:x[0]))
    path = []
    
    for vertice in graph_exon:
        if vertice in ends:
            for end in ends[vertice]:
                append = False
                temp = []
                chief_stack.push(vertice); deputy_stack.push(graph[vertice])
                while chief_stack.stack:
                    adjacency_vertices = deepcopy(deputy_stack.pop())
                    if adjacency_vertices:
                        if len(adjacency_vertices) >= 4:
                            diff_value = [exon[1]-exon[0] for exon in adjacency_vertices]
                            min_index = diff_value.index(min(diff_value))
                            target = adjacency_vertices[min_index]
                            adjacency_vertices = []; adjacency_vertices.append(target)
                        if adjacency_vertices[0] in alternative_edge and chief_stack.get_top() in alternative_front and adjacency_vertices[0] <= end:
                            #判断要添加的exon 是否和前一个 exon 组成可变剪接边
                            if append == False:#此时前面没有可变剪接边
                                temp.append(adjacency_vertices[0])
                                append = True
                            else:#append == True:
                                # if chief_stack.get_top() < temp[-1]:
                                #     continue
                                # else:
                                #此时前面已经有一个可变剪接边 需要判断此时是否能添加另一个可变剪接边
                                if graph_exon.index(adjacency_vertices[0]) - graph_exon.index(temp[-1]) >= 6:
                                    #此时两个可变剪接边的距离大于等于10 则后一个可变剪接边不可以添加
                                    adjacency_vertices.pop(0)
                                elif graph_exon.index(adjacency_vertices[0]) - graph_exon.index(temp[-1]) <= 5:
                                    #此时两个可变剪接边的距离小于等于3 则可以同时存在
                                    if adjacency_vertices[0] not in temp:
                                        temp.append(adjacency_vertices[0])
                        elif adjacency_vertices[0] > end:#如果要添加的可变剪接大于end 则剔除
                            adjacency_vertices.pop(0)
                            
                        if adjacency_vertices:
                            chief_stack.push(adjacency_vertices.pop(0)); deputy_stack.push(adjacency_vertices)
                            deputy_stack.push(graph[chief_stack.get_top()])
                        else:
                            deputy_stack.push(adjacency_vertices)
                            if chief_stack.get_top() == end:
                                path.append(deepcopy(chief_stack.stack))
                                if chief_stack.get_top() in temp:
                                    temp.remove(chief_stack.get_top())
                                    if not temp:
                                        append = False
                                chief_stack.pop(); deputy_stack.pop()
                            else:
                                if chief_stack.get_top() in temp:
                                    temp.remove(chief_stack.get_top())
                                    if not temp:
                                        append = False
                                chief_stack.pop(); deputy_stack.pop()
                    else:
                        deputy_stack.push(adjacency_vertices)
                        if chief_stack.get_top() == end:
                            path.append(deepcopy(chief_stack.stack))
                            if chief_stack.get_top() in temp:
                                temp.remove(chief_stack.get_top())
                                if not temp:
                                    append = False
                            chief_stack.pop(); deputy_stack.pop()
                        else:
                            if chief_stack.get_top() in temp:
                                temp.remove(chief_stack.get_top())
                                if not temp:
                                    append = False
                            chief_stack.pop(); deputy_stack.pop()
    
    return path    


def path_search(splice_graph):
    path = []
    ends = splice_graph.end_point
    
    graph = {}
    for i in splice_graph.vers_inf:
        graph[i] = splice_graph.vers_class[splice_graph.vers_inf.index(i)].adjacency
    
    chief_stack = deepcopy(Stack())
    deputy_stack = deepcopy(Stack())
    
    for vertice in list(graph.keys()):
        if vertice in ends:
            for end in ends[vertice]:
                chief_stack.push(vertice); deputy_stack.push(graph[vertice])
                while chief_stack.stack:
                    adjacency_vertices = deepcopy(deputy_stack.pop())
                    if adjacency_vertices:
                        chief_stack.push(adjacency_vertices.pop(0)); deputy_stack.push(adjacency_vertices)
                        deputy_stack.push(graph[chief_stack.get_top()])
                    else:
                        deputy_stack.push(adjacency_vertices)
                        if chief_stack.get_top() == end:
                            # print(chief_stack.get_top());print(vertice)
                            path.append(deepcopy(chief_stack.stack))
                            chief_stack.pop(); deputy_stack.pop()
                        else:
                            chief_stack.pop(); deputy_stack.pop()
    # print(path)
    return path

def specific_path_search(splice_graph, pair_end):
    path = []
    # pair_end = splice_graph.end_point
    graph = {}
    for i in splice_graph.vers_inf:
        graph[i] = splice_graph.vers_class[splice_graph.vers_inf.index(i)].adjacency
    
    chief_stack = deepcopy(Stack())
    deputy_stack = deepcopy(Stack())
    
    for vertice in list(graph.keys()):
        if vertice in pair_end:
            for end in pair_end[vertice]:
                chief_stack.push(vertice); deputy_stack.push(graph[vertice])
                while chief_stack.stack:
                    adjacency_vertices = deepcopy(deputy_stack.pop())
                    if adjacency_vertices:
                        chief_stack.push(adjacency_vertices.pop(0)); deputy_stack.push(adjacency_vertices)
                        deputy_stack.push(graph[chief_stack.get_top()])
                    else:
                        deputy_stack.push(adjacency_vertices)
                        if chief_stack.get_top() == end:
                            # print(chief_stack.get_top());print(vertice)
                            path.append(deepcopy(chief_stack.stack))
                            chief_stack.pop(); deputy_stack.pop()
                        else:
                            chief_stack.pop(); deputy_stack.pop()
    # print(path)
    return path

# path=path_search(g)

    
if __name__ == '__main__':
    from copy import deepcopy
    #************DFS_test************

    def dfs(vertice,end,graph,trace):
        trace = deepcopy(trace)# 深拷贝，对不同起点，走过的路径不同       
        
        if vertice not in trace:
            trace.append(vertice)
            if vertice == end:
                path.append(trace)
                return None
        
        # 深度优先递归遍历
        for successor_node in graph[vertice]:
            dfs(successor_node,end,graph,trace)

    path = []
    graph = {1: [2,4], 2: [3,4], 3: [4,5], 4: [5], 5:[]}
    ends = {1:[5,4], 3:[5]}

    for k in list(graph.keys()):
        if k in ends:
            for end in ends[k]:
                dfs(k,end,graph,[])
    print(path)
    #************DFS_test************


    #************双栈算法_test*************

    class Stack:
        def __init__(self):
            self.stack = []
            
        def push(self, element):
            self.stack.append(element)
            
        def pop(self):
            return self.stack.pop()
                
        def get_top(self):
            if len(self.stack) > 0:
                return self.stack[-1]
            else:
                return None
            
    chief_stack = deepcopy(Stack())
    deputy_stack = deepcopy(Stack())
    path = []
    graph = {1: [2,4], 2: [3,4], 3: [4,5], 4: [5], 5:[]}
    ends = {1:[5,4], 3:[5]}

    for vertice in list(graph.keys()):
        if vertice in ends:
            for end in ends[vertice]:
                chief_stack.push(vertice); deputy_stack.push(graph[vertice])
                while chief_stack.stack:
                    adjacency_vertices = deepcopy(deputy_stack.pop())
                    if adjacency_vertices:
                        chief_stack.push(adjacency_vertices.pop(0)); deputy_stack.push(adjacency_vertices)
                        deputy_stack.push(graph[chief_stack.get_top()])
                    else:
                        deputy_stack.push(adjacency_vertices)
                        if chief_stack.get_top() == end:
                            # print(chief_stack.get_top());print(vertice)
                            path.append(deepcopy(chief_stack.stack))
                            chief_stack.pop(); deputy_stack.pop()
                        else:
                            chief_stack.pop(); deputy_stack.pop()
    print(path)       
    #************双栈算法_test*************