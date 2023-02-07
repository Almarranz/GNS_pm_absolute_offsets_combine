"""
compare_lists.py
Almarranz 2021
Compare coordinates of two lists and extract the rows whose coordinates columns distance are smaller than variable 'distancia'.
If more than two points in the second list are closer than 'distancia' to a point in the first list,the closest one is choosen 
Input parametres: list1, list2, col_index x1, col_index y1, col_index x2, col_index y2, max distance for a positive match
It computes euclidean distances
Both lists must be numpy arrays
"""
from scipy.spatial import distance
import numpy as np
def compare_lists(list1,list2,ls1_x,ls1_y,ls2_x,ls2_y,distancia):
    ls1_com = np.empty((0,list1.shape[1]))
    ls2_com = np.empty((0,list2.shape[1]))
    coord1 =np.array([list1[:,ls1_x],list1[:,ls1_y]]).T
    coord2 =np.array([list2[:,ls2_x],list2[:,ls2_y]]).T
    for i in range(list1.shape[0]): 
                dist=distance.cdist(coord1[i:i+1,0:2],coord2[:,0:2], 'euclidean')
                d=np.where(dist<distancia)
                if len(d[1])>0:
                    ls1_com = np.append(ls1_com,[list1[i],],axis =0)
                    ls2_com = np.append(ls2_com,[list2[d[1][np.argmin(dist[d])]],],axis =0)
    #dic_listas['stars_1'+str(l)]=diff
    print('common elements in lists ----> ',len(ls1_com))
    return ls1_com,ls2_com