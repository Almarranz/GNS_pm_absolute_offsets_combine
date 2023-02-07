#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 17:06:28 2022

@author: amartinez
"""
import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io.fits import getheader
from astropy.io import fits
from astropy.time import Time
from scipy.spatial import distance
import pandas as pd
import sys
import time
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.spatial import distance
from skimage.transform import estimate_transform
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KDTree
from sklearn.preprocessing import StandardScaler
import os
import math
from scipy.stats import gaussian_kde
import shutil
import astropy.coordinates as ap_coor
from matplotlib.ticker import FormatStrFormatter
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import glob
from datetime import datetime
# %% 
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(10,5)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# 
# %%

pruebas = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_off_comb/pruebas/'
pm_folder ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_off_comb/'
# pm_folder ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_wcs/'

field_one = 7
chip_one = 4
field_two = 7
chip_two =1

align_degree =3#TODO
dmax = 200#TODO
max_sig = 10
vel_cut = 1.4#TODO

# x_dis,y_dis,dvx,dvy,x1,y1,x2,y2,H1,dH1,Ks1,dKs1,H2,dH2,RaH1,DecH1,raH2,decH2
pm_all = np.loadtxt(pm_folder+ 'pm_GaiaRF_ep1_f%sc%s_ep2_f%sc%sdeg%s_dmax%s_sxy%s.txt'%(field_one, chip_one, field_two, chip_two,align_degree,dmax,max_sig))


unc_cut = np.where((pm_all[:,2]<vel_cut) &(pm_all[:,3]<vel_cut))
pm = pm_all[unc_cut]
# mag_cut =np.where((pm[:,12]<19) & (pm[:,12]>11.5))
# pm=pm[mag_cut]
color_cut =np.where((pm[:,8]- pm[:,10])>1.3)
pm=pm[color_cut]

sm=0.1
# sel_m=np.where(abs(pm[:,8]-pm[:,12])<sm)
# pm = pm[sel_m]
# %
bins_ =10
bins_ =30
fig, ax = plt.subplots(1,2,figsize =(20,10))

ax[0].set_title('Gaia RF, Degree = %s'%(align_degree))
ax[1].set_title('Alignment degree = %s, Dmax =%s '%(align_degree,dmax))
# ax[1].set_title('max pm_error =%s , max mag_res =%s, sxy =%s'%(vel_cut,sm,max_sig))

ax[0].hist(pm[:,0], bins =bins_, label = '$\overline{\mu_{l}}$ =%.2f\n$\sigma_{l}$ =%.2f'%(np.mean(pm[:,0]),np.std(pm[:,0])))
ax[0].legend()
ax[1].hist(pm[:,1], bins =bins_, label = '$\overline{\mu_{b}}$ =%.2f\n$\sigma_{b}$ =%.2f'%(np.mean(pm[:,1]),np.std(pm[:,1])))
ax[1].legend()

ax[0].axvline(np.mean(pm[:,0]), color ='red')
ax[1].axvline(np.mean(pm[:,1]), color ='red')
# ax[0].set_xlim(-15,15)
# ax[1].set_xlim(-15,15)
ax[0].set_xlabel('$\mu_{ra}$ (mas/yr)')
ax[1].set_xlabel('$\mu_{dec}$ (mas/yr)')
ax[0].invert_xaxis()

sys.exit('123')
# %%
# DBSCAN
# ra, dec, l, b, pml, pmb,J, H, Ks,x, y, Aks_mean,dAks_mean, radio("), cluster_ID
clust_cand = np.loadtxt('/Users/amartinez/Desktop/PhD/Libralato_data/good_clusters/Sec_B_dmu1_at_minimun_2022-08-30_cluster_num32_2_knn10_area7.49/cluster32_0_0_knn_10_area_7.49_all_color.txt')
cand_coord = SkyCoord(ra = clust_cand[:,0]*u.deg, dec = clust_cand[:,1]*u.deg, frame ='fk5', equinox = 'J2000', obstime = 'J2015.5')

# fig, ax = plt.subplots(1,1)
# ax.scatter(pm[:,14],pm[:,15], color = 'k', alpha= 0.1)
# ax.scatter(clust_cand[:,0],clust_cand[:,1],color = 'lime')
# color = pd.read_csv('/Users/amartinez/Desktop/PhD/python/colors_html.csv')
# strin= color.values.tolist()
# indices = np.arange(0,len(strin),1)


# =============================================================================
# fig,ax = plt.subplots(1,2, figsize =(20,10))
# ax[0].scatter(pm[:,14],pm[:,15])
# =============================================================================
arches_gone = 'no'#TODO 'yes' of you want to remove the arches from the data. If not it only finds the arches.
if arches_gone == 'yes':
    #   0      1    2  3  4   5  6  7 8   9   10  11  12  13  14    15   16    17
    # x_dis,y_dis,dvx,dvy,x1,y1,x2,y2,H1,dH1,Ks1,dKs1,H2,dH2,RaH1,DecH1,raH2,decH2
    # ar_g = np.where((pm[:,14]>266.45) & (pm[:,14]<266.47)
    #                   & (pm[:,15]>-28.83) & (pm[:,15]<-28.82))
    m_arc = SkyCoord(ra = [266.46107138], dec = [-28.82241165],unit = 'degree' )
    all_coor =  SkyCoord(ra = pm[:,-2], dec = pm[:,-1],unit = 'degree' )
    idxc, group_md, d2d,d3d =  ap_coor.search_around_sky(m_arc, all_coor, 20*u.arcsec)
    
    pm =np.delete(pm,group_md,axis=0)
# =============================================================================
# ax[1].scatter(pm[:,14],pm[:,15])   
# =============================================================================
    
divis_list = [1]#TODO
# divis_list = [1]#TODO

# samples_list =[10,9,8,7,6,5]#TODO
samples_list =[20]#TODO
sim_lim = 'mean'#TODO options: minimun or mean
gen_sim ='Shuffle'#TODO it is not yet implemented shuffle
clustered_by = 'pm_color'#TODO
# clustered_by = 'pm'#TODO
clus_num = 0
for div in range(len(divis_list)):
    divis = divis_list[div]
    xg = np.linspace(min(pm[:,4]),max(pm[:,4]),(divis+1)*2-1)
    yg = np.linspace(min(pm[:,5]),max(pm[:,5]),(divis+1)*2-1)
    for sd in range(len(samples_list)):
        samples_dist = samples_list[sd]
        for xi in range(len(xg)-2):
            for yi in range(len(xg)-2):
        # for xi in range(3):
        #     for yi in range(1):
                # fig, ax = plt.subplots(1,1)
                # ax.scatter(pm[:,4],pm[:,5],color ='k',alpha=0.1)
                valid = np.where((pm[:,4] < xg[xi+2])
                                 & (pm[:,4] >xg[xi]) 
                                 & (pm[:,5] < yg[yi+2]) 
                                 & (pm[:,5] >yg[yi]))
                pm_sub = pm[valid]
                # ax.scatter(pm_sub[:,6],pm_sub[:,7],s =5,color=strin[np.random.choice(indices)] )
                coordenadas = SkyCoord(ra=pm_sub[:,14]*u.degree, dec=pm_sub[:,15]*u.degree,frame ='icrs', equinox = 'J2000', obstime = 'J2022.4')#
                mul,mub = pm_sub[:,0],pm_sub[:,1]
                x,y = pm_sub[:,6], pm_sub[:,7]
                colorines = pm_sub[:,8]-pm_sub[:,10]
                H_datos, K_datos = pm_sub[:,8], pm_sub[:,10]
                
                area =  (xg[xi+1]- xg[xi])*(yg[yi+1]- yg[yi])/1000**2/3600
                mul_kernel, mub_kernel = gaussian_kde(mul), gaussian_kde(mub)
                x_kernel, y_kernel = gaussian_kde(x), gaussian_kde(y)
                color_kernel = gaussian_kde(colorines)
        
                if clustered_by =='pm_color':
                   X=np.array([mul,mub,x,y,colorines]).T
                   X_stad = StandardScaler().fit_transform(X)
                   tree = KDTree(X_stad, leaf_size=2) 
                if clustered_by == 'pm':
                   X=np.array([mul,mub,x,y]).T
                   X_stad = StandardScaler().fit_transform(X)
                   tree = KDTree(X_stad, leaf_size=2) 
               
                dist, ind = tree.query(X_stad, k=samples_dist) #DistNnce to the 1,2,3...k neighbour
                d_KNN=sorted(dist[:,-1])#distance to the Kth neighbour
                lst_d_KNN_sim = []
                for d in range(5):
                    if clustered_by =='pm_color':
                        mub_sim,  mul_sim = mub_kernel.resample(len(pm_sub)), mul_kernel.resample(len(pm_sub))
                        x_sim, y_sim = x_kernel.resample(len(pm_sub)), y_kernel.resample(len(pm_sub))
                        color_sim = color_kernel.resample(len(pm_sub))
                        X_sim=np.array([mul_sim[0],mub_sim[0],x_sim[0],y_sim[0],color_sim[0]]).T
                        X_stad_sim = StandardScaler().fit_transform(X_sim)
                        tree_sim =  KDTree(X_stad_sim, leaf_size=2)
                        
                        dist_sim, ind_sim = tree_sim.query(X_stad_sim, k=samples_dist) #DistNnce to the 1,2,3...k neighbour
                        d_KNN_sim=sorted(dist_sim[:,-1])#distance to the Kth neighbour
                    
                        lst_d_KNN_sim.append(min(d_KNN_sim))    
                    if clustered_by == 'pm':
                        mub_sim,  mul_sim = mub_kernel.resample(len(pm_sub)), mul_kernel.resample(len(pm_sub))
                        x_sim, y_sim = x_kernel.resample(len(pm_sub)), y_kernel.resample(len(pm_sub))
                        X_sim=np.array([mul_sim[0],mub_sim[0],x_sim[0],y_sim[0]]).T
                        X_stad_sim = StandardScaler().fit_transform(X_sim)
                        tree_sim =  KDTree(X_stad_sim, leaf_size=2)
                        
                        dist_sim, ind_sim = tree_sim.query(X_stad_sim, k=samples_dist) #DistNnce to the 1,2,3...k neighbour
                        d_KNN_sim=sorted(dist_sim[:,-1])#distance to the Kth neighbour
                        
                        lst_d_KNN_sim.append(min(d_KNN_sim))
                d_KNN_sim_av = np.mean(lst_d_KNN_sim)
# =============================================================================
#                 fig, ax = plt.subplots(1,1,figsize=(10,10))
#                 ax.hist(d_KNN,bins ='auto',histtype ='step',color = 'k')
#                 ax.hist(d_KNN_sim,bins ='auto',histtype ='step',color = 'r')
#                 ax.set_xlabel('%s-NN distance'%(samples_dist)) 
#                 
# =============================================================================
                if sim_lim == 'mean':
                    eps_av = round((min(d_KNN)+d_KNN_sim_av)/2,3)
                    valor = d_KNN_sim_av
                elif sim_lim == 'minimun':
                    eps_av = round((min(d_KNN)+min(lst_d_KNN_sim))/2,3)
                    valor = min(lst_d_KNN_sim)
                texto = '\n'.join(('min real d_KNN = %s'%(round(min(d_KNN),3)),
                                'min sim d_KNN =%s'%(round(valor,3)),
                                'average = %s'%(eps_av),'%s'%(sim_lim),'%s'%(gen_sim)))
                
# =============================================================================
#                 props = dict(boxstyle='round', facecolor='w', alpha=0.5)
#                 # place a text box in upper left in axes coords
#                 ax.text(0.55, 0.25, texto, transform=ax.transAxes, fontsize=20,
#                     verticalalignment='top', bbox=props)
#                 
#                 ax.set_ylabel('N') 
# =============================================================================
                
               
                clustering = DBSCAN(eps=eps_av, min_samples=samples_dist).fit(X_stad)
                l=clustering.labels_
                
                n_clusters = len(set(l)) - (1 if -1 in l else 0)
                # print('Group %s.Number of cluster, eps=%s and min_sambles=%s: %s'%(group,round(epsilon,2),samples,n_clusters))
                n_noise=list(l).count(-1)
                # %
                u_labels = set(l)
                colors=[plt.cm.rainbow(i) for i in np.linspace(0,1,len(set(l)))]# Returns a color for each cluster. Each color consists in four number, RGBA, red, green, blue and alpha. Full opacity black would be then 0,0,0,1
                # %
                
                # %
                for k in range(len(colors)): #give noise color black with opacity 0.1
                    if list(u_labels)[k] == -1:
                        colors[k]=[0,0,0,0.1]
                # %      
                colores_index=[]
                
                for c in u_labels:
                    cl_color=np.where(l==c)
                    colores_index.append(cl_color)
                
                for i in range(len(set(l))-1):
                    # ax[0].scatter(pm_sub[:,2][colores_index[i]], pm_sub[:,3][colores_index[i]], color = 'orange')
                    c2 = SkyCoord(ra = pm_sub[:,16][colores_index[i]],dec = pm_sub[:,17][colores_index[i]], unit ='degree',  equinox = 'J2000', obstime = 'J2015.4')
                    fig, ax = plt.subplots(1,3,figsize=(30,10))
                    color_de_cluster = 'lime'
                    # fig, ax = plt.subplots(1,3,figsize=(30,10))
                    # ax[2].invert_yaxis()
                   
                    ax[0].set_title('Min %s-NN= %s. cluster by: %s '%(samples_dist,round(min(d_KNN),3),clustered_by))
                    # t_gal['l'] = t_gal['l'].wrap_at('180d')
                    ax[0].scatter(X[:,0][colores_index[-1]],X[:,1][colores_index[-1]], color=colors[-1],s=50,zorder=1)
                    ax[0].scatter(X[:,0],X[:,1], color=colors[-1],s=50,zorder=1)
                    # ax[1].quiver(t_gal['l'][colores_index[-1]].value,t_gal['b'][colores_index[-1]].value, X[:,0][colores_index[-1]]-pms[2], X[:,1][colores_index[-1]]-pms[3], alpha=0.5, color=colors[-1])
            
                    ax[0].scatter(X[:,0][colores_index[i]],X[:,1][colores_index[i]], color=color_de_cluster ,s=50,zorder=3)
                    # ax[0].set_xlim(-10,10)
                    # ax[0].set_ylim(-10,10)
                    ax[0].set_xlabel(r'$\mathrm{\mu_{ra} (mas\ yr^{-1})}$',fontsize =30) 
                    ax[0].set_ylabel(r'$\mathrm{\mu_{dec} (mas\ yr^{-1})}$',fontsize =30) 
                    ax[0].invert_xaxis()
                    
                    mul_sig, mub_sig = np.std(X[:,0][colores_index[i]]), np.std(X[:,1][colores_index[i]])
                    mul_mean, mub_mean = np.mean(X[:,0][colores_index[i]]), np.mean(X[:,1][colores_index[i]])
                    
                    mul_sig_all, mub_sig_all = np.std(X[:,0]), np.std(X[:,1])
                    mul_mean_all, mub_mean_all = np.mean(X[:,0]), np.mean(X[:,1])
                     
                    ax[0].axvline(mul_mean_all, color ='red')
                    ax[0].axhline(mub_mean_all, color ='red')
                    
                    ax[0].axvline(mul_mean_all + mul_sig_all,linestyle = 'dashed', color ='red')
                    ax[0].axvline(mul_mean_all - mul_sig_all,linestyle = 'dashed', color ='red')
                    
                    ax[0].axhline(mub_mean_all + mub_sig_all,linestyle = 'dashed', color ='red')
                    ax[0].axhline(mub_mean_all - mub_sig_all,linestyle = 'dashed', color ='red')
                
                    vel_txt = '\n'.join((r'$\mathrm{\mu_{ra}}$ = %s, $\mathrm{\mu_{dec}}$  = %s'%(round(mul_mean,3), round(mub_mean,3)),
                                         '$\sigma_{ra}$ = %s, $\sigma_{dec}$ = %s'%(round(mul_sig,3), round(mub_sig,3)))) 
                    vel_txt_all = '\n'.join(('mura = %s, mudec = %s'%(round(mul_mean_all,3), round(mub_mean_all,3)),
                                         '$\sigma_{ra}$ = %s, $\sigma_{dec}$ = %s'%(round(mul_sig_all,3), round(mub_sig_all,3))))
                    
                    propiedades = dict(boxstyle='round', facecolor=color_de_cluster , alpha=0.2)
                    propiedades_all = dict(boxstyle='round', facecolor=colors[-1], alpha=0.1)
                    ax[0].text(0.05, 0.95, vel_txt, transform=ax[0].transAxes, fontsize=30,
                        verticalalignment='top', bbox=propiedades)
                    ax[0].text(0.05, 0.15, vel_txt_all, transform=ax[0].transAxes, fontsize=20,
                        verticalalignment='top', bbox=propiedades_all)   
                    #This calcualte the maximun distance between cluster members to have a stimation of the cluster radio
                    sep = [max(c2[c_mem].separation(c2)) for c_mem in range(len(c2))]
                    rad = max(sep)/2
                    
                    m_point = SkyCoord(ra =[np.mean(c2.ra)], dec = [np.mean(c2.dec)],frame ='icrs', equinox = 'J2000', obstime = 'J2022.4')             
                    print(m_point)
                    
                    idxc, group_md, d2d,d3d =  ap_coor.search_around_sky(m_point,coordenadas, rad*2)
                   
                    ax[0].scatter(mul[group_md],mub[group_md], color='red',s=50,zorder=1,marker='x',alpha = 0.7)
        
                    prop = dict(boxstyle='round', facecolor=color_de_cluster , alpha=0.2)
                    ax[1].text(0.15, 0.15, 'aprox cluster radio = %s"\n cluster stars = %s '%(round(rad.to(u.arcsec).value,2),len(colores_index[i][0])), transform=ax[1].transAxes, fontsize=30,
                                            verticalalignment='top', bbox=prop)
                    # ax[1].scatter(MS_coord.ra, MS_coord.dec, s=20, color ='b', marker ='.')
                    ax[1].scatter(pm_sub[:,14], pm_sub[:,15], color='k',s=50,zorder=1,alpha=0.1)#
                    # ax[1].scatter(datos[:,5],datos[:,6],color='k' ,s=50,zorder=1,alpha=0.01)
                    ax[1].scatter(pm_sub[:,14][colores_index[i]],pm_sub[:,15][colores_index[i]],color=color_de_cluster ,s=50,zorder=3)
        
                    ax[1].scatter(pm_sub[:,14][group_md],pm_sub[:,15][group_md],s=50,color='r',alpha =0.1,marker ='x')
                    ax[1].set_xlabel('Ra(deg)',fontsize =30) 
                    ax[1].set_ylabel('Dec(deg)',fontsize =30) 
                    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    # ax[1].set_title('Area_%s_%s =%.2f $arcmin^{2}$'%(xi,yi,area))
                    ax[1].set_title('Degree = %s'%(align_degree))
        
                    
                    for cand_star in range(len(cand_coord)):
                        sep_cand = c2.separation(cand_coord[cand_star])
                        if min(sep_cand.value) <1/3600:
                            print('SOMETHING is CLOSE!')
                            ax[1].scatter(cand_coord.ra, cand_coord.dec,color ='b',s =80,marker ='x',zorder =3)
                            ax[0].set_facecolor('lavender')
                            ax[1].set_facecolor('lavender')
                            ax[2].set_facecolor('lavender')
                    
                    ax[2].set_title('clus_num = %s,i =%s,knn = %s,area =%.2f'%(clus_num,i,samples_dist,area))
                    ax[2].scatter(pm_sub[:,8]-pm_sub[:,10],pm_sub[:,10],alpha=0.1,color ='k')  
                    ax[2].scatter(pm_sub[:,8][group_md]-pm_sub[:,10][group_md],pm_sub[:,10][group_md],alpha=0.7,c='r',marker = 'x')
                    txt_around = '\n'.join(('H-Ks =%.3f'%(np.median(pm_sub[:,8][group_md]-pm_sub[:,10][group_md])),
                                         '$\sigma_{H-Ks}$ = %.3f'%(np.std(pm_sub[:,8][group_md]-pm_sub[:,10][group_md])),
                                         'diff_color = %.3f'%(max(pm_sub[:,8][group_md]-pm_sub[:,10][group_md])-min(pm_sub[:,8][group_md]-pm_sub[:,10][group_md]))))
                    props_arou = dict(boxstyle='round', facecolor='r', alpha=0.3)
                    ax[2].text(0.50, 0.25,txt_around, transform=ax[2].transAxes, fontsize=30,
                        verticalalignment='top', bbox=props_arou)       
                    ax[2].scatter(pm_sub[:,8][colores_index[i]]-pm_sub[:,10][colores_index[i]],pm[:,10][colores_index[i]],alpha=1,color =color_de_cluster,s=100)
                    ax[2].invert_yaxis()  
                    ax[2].set_xlabel('$H-Ks$',fontsize =30)
                    ax[2].set_ylabel('$Ks$',fontsize =30)
                    ax[2].set_xlim(1,3)
                    txt_color = '\n'.join(('H-Ks =%.3f'%(np.median(pm_sub[:,8][colores_index[i]]-pm_sub[:,10][colores_index[i]])),
                                                            '$\sigma_{H-Ks}$ = %.3f'%(np.std(pm_sub[:,8][colores_index[i]]-pm_sub[:,10][colores_index[i]])),
                                                            'diff_color = %.3f'%(max(pm_sub[:,8][colores_index[i]]-pm_sub[:,10][colores_index[i]])-min(pm_sub[:,8][colores_index[i]]-pm_sub[:,10][colores_index[i]]))))
                    props = dict(boxstyle='round', facecolor=color_de_cluster, alpha=0.3)
                    ax[2].text(0.50, 0.95, txt_color, transform=ax[2].transAxes, fontsize=30,
                                            verticalalignment='top', bbox=props)
                    plt.show()
                    
# =============================================================================
#                     emilio = '/Users/amartinez/Desktop/PhD/for_Emilio/'
#                     save_for_emilio =input('save for emilio?')
#                     if save_for_emilio == 'yes':
#                         print(clus_num,i)
#                         np.savetxt(emilio + 'clusters_gns_arches_%s_%s.txt'%(clus_num,i),np.c_[pm_sub[:,14][group_md],pm_sub[:,15][group_md],X[group_md]], fmt ='%.8f',
#                                    header ='ra, dec, pmra, pmdec, x, y, color')
# =============================================================================
                        
                    clus_array = np.array([pm_sub[:,14][colores_index[i]],pm_sub[:,15][colores_index[i]],
                                          X[:,0][colores_index[i]],X[:,1][colores_index[i]],
                                          X[:,2][colores_index[i]],X[:,3][colores_index[i]],
                                          H_datos[colores_index[i]], K_datos[colores_index[i]]]).T
                    # =============================================================================
                            #             Here it compare the cluster you want to save wiith the rest of the 
                            #             saved cluster if repited, it saves in the same cluster 
                            #             
                            # =============================================================================
                                         
                    frase = 'Do you want to save this cluster?'
                    print('\n'.join((len(frase)*'π',frase+'\n("yes" or "no")',len(frase)*'π')))
                    save_clus = input('Awnser:')
                    # save_clus = 'yes'
                    # np.savetxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/Sec_7_clus/Arches_gns_gRF/Arches_gns_gRF.txt',clus_array,fmt='%.8f ',header ='ra, dec, pml, pmb ,x,y, H, Ks')
                    print('You said: %s'%(save_clus))
                    if save_clus =='yes' or save_clus =='y':
                        
                        intersection_lst =[]
                        section = field_two
                        ic =xi
                        jr =yi
                        check_folder = glob.glob(pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one)+'cluster_num*')
                        if len(check_folder) == 0:
                            os.makedirs(pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one) +'cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area))
                            np.savetxt(pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one)+'cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area)+
                                       'cluster%s_%.0f_%.0f_knn_%s_area_%.2f_%s.txt'%(clus_num,xi,yi,samples_dist,area,clustered_by),clus_array,
                                       fmt='%.8f ',
                                       header ='ra, dec, pml, pmb ,x,y, H, Ks')
                            ax[2].set_title('Saved in cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area))
                            # plt.show()
                            clus_num +=1   
                        else:
                            break_out_flag = False
                            for f_check in check_folder:
                                clus_lst = os.listdir(f_check)
                                for n_txt in clus_lst:
                                    ra_dec = np.loadtxt(f_check+'/'+n_txt,usecols=(0,1))
                                    ra_dec_clus = clus_array[:,0:2]
                                    aset = set([tuple(x) for x in ra_dec_clus])
                                    bset = set([tuple(x) for x in ra_dec])
                                    intersection = np.array([x for x in aset & bset])
                                    intersection_lst.append(len(intersection))
                                    # print('This is intersection',intersection_lst)
                                    if len(intersection)> 0 :
                                        print('Same (or similar) cluster  is in %s'%(f_check))
                                        np.savetxt(f_check+'/'+
                                                   'cluster%s_%.0f_%.0f_knn_%s_area_%.2f_%s.txt'%(clus_num,ic,jr,samples_dist,area,clustered_by),clus_array,
                                                   fmt='%.8f ',
                                                   header ='ra, dec, pml, pmb ,x,y, H, Ks')
                                        ax[2].set_title('Saved in %s'%(os.path.basename(f_check)))
                                        # plt.show()
                                        clus_num +=1 
                                        break_out_flag = True
                                        break
                                if break_out_flag:
                                    break
                                
                            if np.all(np.array(intersection_lst)==0):
                                # clus_num +=1
                                print('NEW CLUSTER')
                                os.makedirs(pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one) +'cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area))
                                np.savetxt(pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one) +'cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area)+
                                           'cluster%s_%.0f_%.0f_knn_%s_area_%.2f_%s.txt'%(clus_num,xi/0.5,yi/0.5,samples_dist,area,clustered_by),clus_array,
                                           fmt='%.8f ',
                                           header ='ra, dec, pml, pmb ,x,y, H, Ks')
                                ax[2].set_title('Saved in cluster_num%s_%s_knn%s_area%.2f/'%(clus_num,i,samples_dist,area))
                                # plt.show()
                                clus_num +=1   
                                   
                                # read_txt = glob.glob(check_folder[f_check]+'/cluster_*')
                                # for clust_text in range(len(read_txt)):
                                #     print(read_txt[clust_text])
                                
                      
                    elif save_clus =='stop':
                        frase = 'Do you want to copy the folder with the clusters into the pruebas directory?\n("yes" or "no")'
                        print('\n'.join((len(frase)*'',frase,len(frase)*'')))
                        save_folder = input('Awnser:')   
                        if save_folder == 'yes' or save_folder == 'y':       
                            source_dir = pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one)
                            # destination_dir = '/Users/amartinez/Desktop/morralla/Sec_%s_at_%s'%(section,datetime.now())
                            destination_dir = pruebas +'Sec_f%sc%s_at_%s'%(field_one, chip_one,datetime.now())
                            shutil.copytree(source_dir, destination_dir)
                            sys.exit('You stoped it')
                        else:
                            sys.exit('You stoped it')
                        sys.exit('Chao')
                    
                    else:
                        print('No saved')
# %%               
mix_color = 'yes'
dmu_lim = vel_cut
frase = 'Do you want to copy the folder with the clusters into the morralla directory?\n("yes" or "no")'
print('\n'.join((len(frase)*'',frase,len(frase)*'')))
save_folder = input('Awnser:')   
if save_folder == 'yes' or save_folder == 'y':       
    source_dir = pruebas + 'Sec_f%sc%s_clus/'%(field_one,chip_one)
    if mix_color == 'yes':
        destination_dir = '/Users/amartinez/Desktop/morralla/GNS1_f%sc%s_dmu%s_at_GNS2_f%sc%s_%s'%(field_one,chip_one,vel_cut,field_two,chip_two,datetime.now())
    elif mix_color == 'no':
        destination_dir = '/Users/amartinez/Desktop/morralla/Sec_%s_dmu%s_at_%s_%s_nocolor_%s'%(section,dmu_lim,sim_lim,gen_sim,datetime.now())
    shutil.copytree(source_dir, destination_dir)
    sys.exit('You stoped it')
else:
    sys.exit('You stoped it')
sys.exit('Chiao')
# %%
# divis_list = [5]#TODO


# # samples_list =[10,9,8,7,6,5]#TODO
# samples_list =[5]#TODO
# sim_lim = 'mean'#TODO options: minimun or mean
# gen_sim ='Kernnel'#TODO it is not yet implemented shuffle
# clustered_by = 'pm_color'#Reminiscent of a previous script
# cluster_by = 'pm'
# clus_num = 0
# for div in range(len(divis_list)):
#     divis = divis_list[div]
#     xg = np.linspace(min(pm[:,4]),max(pm[:,4]),(divis+1)*2-1)
#     yg = np.linspace(min(pm[:,5]),max(pm[:,5]),(divis+1)*2-1)
#     for sd in range(len(samples_list)):
#         samples_dist = samples_list[sd]
#         for xi in range(len(xg)-2):
#             for yi in range(len(xg)-2):
#         # for xi in range(3):
#         #     for yi in range(1):
#                 fig, ax = plt.subplots(1,1)
#                 ax.scatter(pm[:,4],pm[:,5],color ='k',alpha=0.1)
#                 valid = np.where((pm[:,4] < xg[xi+2])
#                                  & (pm[:,4] >xg[xi]) 
#                                  & (pm[:,5] < yg[yi+2]) 
#                                  & (pm[:,5] >yg[yi]))
#                 pm_sub = pm[valid]
#                 ax.scatter(pm_sub[:,6],pm_sub[:,7],s =5,color=strin[np.random.choice(indices)] )








