from laspy.file import File
import numpy as np
from scipy import spatial
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import pgSQL_3DCityDB_reader as pgSQL
import time
import math
from shapely.geometry import Polygon, Point

def data_loader(input_file):
    print(input_file)
    inFile = File(input_file, mode= 'rw')
    for ind in range(0,len(inFile.x)):
        inFile.user_data[ind] = 0
    #print(inFile.X,inFile.Y,inFile.Z)
    #print(inFile.x, inFile.y, inFile.z)
    #tree = spatial.cKDTree(list(zip(inFile.x, inFile.y, inFile.z)))
    return inFile, []

def crop_pc(inFile, bbox, id):
    try:
        outFile = File("cropped_buildings/building_id_{}.las".format(id), mode= 'rw')
    except:
        x_invalid = np.logical_or((bbox[0][0]-1 >= inFile.x),(bbox[1][0]+1 <= inFile.x))
        y_invalid = np.logical_or((bbox[0][1]-1 >= inFile.y),(bbox[1][1]+1 <= inFile.y))
        #z_invalid = np.logical_or((bbox[0][2] >= inFile.z),(bbox[1][2] <= inFile.z))
        indices = np.logical_not(np.logical_or(x_invalid, y_invalid))
        #print(True in indices)
        if True not in indices:
            print("empty case")
            return "none", "none"
        outFile = File("cropped_buildings/building_id_{}.las".format(id), mode='w', header=inFile.header)
        outFile.points = inFile.points[indices]
    tree = spatial.cKDTree(list(zip(outFile.x, outFile.y, outFile.z)))
    return outFile, tree

def data_h_filter(input_file):
    print(input_file)
    inFile = File(input_file, mode='rw')
    filter = inFile.z > 4.5
    outFile = File(input_file.replace(".las",'')+"height_corr.las", mode='w', header=inFile.header)
    outFile.points=inFile.points[filter]
    for ind in range(0,len(outFile.x)):
        outFile.user_data[ind] = 0
    tree = spatial.cKDTree(list(zip(outFile.x, outFile.y, outFile.z)))
    return outFile, tree

def comp_normals(inFile, tree, tolerance):
    #coord_list = list(zip(inFile.X, inFile.Y))
    for index in range(0,len(inFile.X)):
        dist, neighbors = tree.query([inFile.x[index], inFile.y[index],inFile.z[index]], 10)
        neighb_coord = []
        for ind in neighbors:
            neighb_coord.append([inFile.x[ind], inFile.y[ind], inFile.z[ind]])
        normal, ref_pt = fit_plane(neighb_coord)
        #print(normal)
        #if normal[2]<0:
            #print(normal)
            #normal=normal*-1
            #print(normal)
        #inFile.user_data[index] = normal[2]/(abs(normal[1])+abs(normal[0]))
        if abs(normal[2])>(tolerance*abs(normal[1])) and abs(normal[2])>(tolerance*abs(normal[0])):
        #if normal[2] > (tolerance * (normal[1]+normal[0])):
            inFile.user_data[index] = 1
        else:
            inFile.user_data[index] = 0

def comp_height_diff(inFile, tree):
    #coord_list = list(zip(inFile.X, inFile.Y))
    height_dev = []
    for index in range(0,len(inFile.x)):
        dist, neighbors = tree.query([inFile.x[index], inFile.y[index],inFile.z[index]], 3)
        neighb_z = []
        for ind in neighbors:
            neighb_z.append(abs(inFile.z[index] - inFile.z[ind]))
        #print(max(neighb_z))
        if max(neighb_z)>0.05:
            height_dev.append(index)
            inFile.user_data[index] = 1
        else:
            inFile.user_data[index] = 0
    return height_dev

def fit_plane(pts):
    """
    Fits a plane through a set of points using principal component analysis.

    Input:
        pts:    the points to fit a plane through

    Output:
        (n,c)   a tuple with a point p that lies on the plane and the normalised normal vector n of the plane.
    """
    #print(pts)
    # shift points to mean
    mean = np.mean(pts, axis=0)
    pts -= mean
    # compute covariance matrix and eigenvalues and eignevectors
    cov = np.cov(pts, rowvar=False)
    evals, evecs = eigh(cov)
    # find smallest eigenvalue, the corresponging eigenvector is our normal n
    idx = np.argsort(evals)[::-1]
    #print(idx, evals)
    evecs = evecs[:, idx]
    evals = evals[idx]
    #print(evecs)
    n = evecs[:,-1]
    if n[2]<0:
        n = n*-1
    #c = mean
    curvature = evals[0]/(np.sum(evals))
    #print("curvature test", evals[0], curvature)
    #print(n)

    return n, curvature

#deprecated
def height_deviations(inFile, min_dev, extents_min=(90000.6570000001, 431957.705000001, 10.95), extents_max=(90041.332, 432069.438000002, 10.95)):
    dev_points = []
    for index in range(0,len(inFile.x)):
        #if 0.2 < abs(inFile.z[index] - extents_min[2]) < 2.0:
            #print(inFile.x[index], inFile.y[index],inFile.z[index])
        if extents_min[0]<inFile.x[index]<extents_max[0] and extents_min[1]<inFile.y[index]<extents_max[1]:
            #print(inFile.Z - extents_min[2])
            if min_dev<abs(inFile.z[index]-extents_min[2])<2.0:
                dev_points.append(index)
                inFile.user_data[index] = 1
            else:
                #inFile.user_data[index] = 0
                pass
        else:
            #print("point outside")
            #inFile.user_data[index] = 0
            pass
    return dev_points

def height_deviations_normal(inFile, min_dist, pt_poly_dict, surf_list, normals):
    dev_seeds = []
    for surf_idx in range(0,len(surf_list)):
        curr_surf_seeds = []
        curr_ref_pt = surf_list[surf_idx][0]
        for pt_idx in pt_poly_dict[surf_idx]:
            distance = shortest_distance(curr_ref_pt, normals[surf_idx], (inFile.x[pt_idx], inFile.y[pt_idx], inFile.z[pt_idx]))
            if distance > min_dist:
                curr_surf_seeds.append(pt_idx)
                inFile.user_data[pt_idx] = 10
            else:
                inFile.user_data[pt_idx] = 0
        dev_seeds.append(curr_surf_seeds)
    return dev_seeds

def height_deviations_vertical(inFile, min_dist, pt_poly_dict, surf_list, normals):
    dev_seeds = []
    for surf_idx in range(0,len(surf_list)):
        curr_surf_seeds = []
        curr_ref_pt = surf_list[surf_idx][0]
        for pt_idx in pt_poly_dict[surf_idx]:
            shortest_vector_length = shortest_distance(curr_ref_pt, normals[surf_idx], (inFile.x[pt_idx], inFile.y[pt_idx], inFile.z[pt_idx]))
            curr_normal = normals[surf_idx]*shortest_vector_length
            displacement_v = [np.nan, np.nan, np.nan]
            displacement_v[0] = - curr_normal[0]
            displacement_v[1] = - curr_normal[1]
            displacement_v[2] = -(displacement_v[0]*curr_normal[0]+displacement_v[1]*curr_normal[1])/curr_normal[2]
            distance = displacement_v[2]+curr_normal[2]
            #print(distance)
            if distance > min_dist:
                curr_surf_seeds.append(pt_idx)
                inFile.user_data[pt_idx] = 10
            else:
                inFile.user_data[pt_idx] = 0
        #print((normals[surf_idx][0]**2+normals[surf_idx][1]**2+normals[surf_idx][2]**2)**2)
        dev_seeds.append(curr_surf_seeds)
    return dev_seeds

def shortest_distance(plane_pt, normal_v, query_pt):
    plane_query_v = np.array([plane_pt[0]-query_pt[0],plane_pt[1]-query_pt[1],plane_pt[2]-query_pt[2]])
    #print(abs(np.dot(plane_query_v,normal_v)))
    return abs(np.dot(plane_query_v,normal_v))

def map_curvature(inFile, tree, neighb_num):
    for index in range(0,len(inFile.X)):
        dist, neighbors = tree.query([inFile.x[index], inFile.y[index],inFile.z[index]], neighb_num)
        neighb_coord = []
        for ind in neighbors:
            neighb_coord.append([inFile.x[ind], inFile.y[ind], inFile.z[ind]])
        normal, curvature = fit_plane(neighb_coord)
        #print(int(round(curvature*10)), curvature)
        inFile.user_data[index] = int(round(curvature*100))

def neighb_filter(neighb_inds, curr_region, region_dict):
    outlist = []
    for neighb_ind in neighb_inds:
        if (neighb_ind in curr_region) == False and (neighb_ind in region_dict) == False:
            outlist.append(neighb_ind)
    return outlist

def grow_regions(deviations, tree, inFile, tolerance):
    print("num deviations", len(deviations))
    SeedInd = deviations[::10] #select every 10th point from the deviations as a seed
    #for lookup speed, convert list into dictionnary
    deviations_dict = {}
    for ind in deviations:
        deviations_dict[ind] = 1
    region_dict = {}
    seed_dict = {}
    # for ind in SeedInd:
    #     seed_dict[ind] = 1
    reg_idx = 0
    stop_code = False
    add_pts_cnt = 0
    output_regions = []
    while(stop_code==False):
        SeedPtInd = SeedInd[reg_idx]
        # in case the previous iteration arrived at the last seed, search for new seeds
        if reg_idx == len(SeedInd)-1:
            print(reg_idx, len(SeedInd), "looking for new seeds")
            pt_add_counter = 0
            for index in range(0, len(deviations)):
                #print("looking for orphan points")
                if deviations[index] not in region_dict:
                    SeedInd.append(deviations[index])
                    pt_add_counter += 1
            # num_new_seeds = int(len(unclassified))
            # print('%.20f' % len(unclassified))
            if pt_add_counter == 0:
                print("no more seeds, exit this region growing")
                print("number of added points:", add_pts_cnt)
                break
            print(len(SeedInd))
        #print(reg_idx)
        stack = [SeedPtInd]
        curr_region = {}
        #check if the seed point is already classified
        if SeedPtInd in region_dict:
            #print("seed already classified")
            reg_idx +=1
            continue
        #curr_region[SeedPtInd] = 1
        while len(stack)>0:
            exp_pt_ind = stack.pop()
            #print("popped point")
            exp_pt = [inFile.x[exp_pt_ind], inFile.y[exp_pt_ind], inFile.z[exp_pt_ind]]
            #now let's find neighbors
            dist, neighb = tree.query(exp_pt, 10)
            neighb = neighb_filter(neighb, curr_region, region_dict)
            for n_pt_ind in neighb:
                if n_pt_ind in deviations_dict:
                    #print("added point")
                    curr_region[n_pt_ind] = reg_idx+1
                    stack.append(n_pt_ind)
                    #print(n_pt_ind)
                    pass
                else:
                    n_pt = [inFile.x[n_pt_ind], inFile.y[n_pt_ind], inFile.z[n_pt_ind]]
                    dist, plane_neighb = tree.query(n_pt, 10)
                    plane_pts = []
                    for plane_ind in plane_neighb:
                        plane_pts.append([inFile.x[plane_ind], inFile.y[plane_ind], inFile.z[plane_ind]])
                    normal, curvature = fit_plane(plane_pts)
                    # # coord_list = list(zip(inFile.X, inFile.Y))
                    # for index in range(0, len(inFile.X)):
                    #     dist, neighbors = tree.query([inFile.x[index], inFile.y[index], inFile.z[index]], 10)
                    #     neighb_coord = []
                    #     for ind in neighbors:
                    #         neighb_coord.append([inFile.x[ind], inFile.y[ind], inFile.z[ind]])
                    #     normal, ref_pt = fit_plane(neighb_coord)
                    #     # print(normal)
                    #     # if normal[2]<0:
                    #     # print(normal)
                    #     # normal=normal*-1
                    #     # print(normal)
                    #     # inFile.user_data[index] = normal[2]/(abs(normal[1])+abs(normal[0]))
                    #     if abs(normal[2]) > (tolerance * abs(normal[1])) and abs(normal[2]) > (
                    #             tolerance * abs(normal[0])):
                    #         # if normal[2] > (tolerance * (normal[1]+normal[0])):
                    #         inFile.user_data[index] = 1
                    #     else:
                    #         inFile.user_data[index] = 0
                    if abs(normal[2]) > (tolerance * abs(normal[1])) and abs(normal[2]) > (tolerance * abs(normal[0])):
                        #curr_region[n_pt_ind] = -1
                        #print("rejected point after plane check")
                        pass
                    else:
                        #print(normal, curvature, "added point after normal vector check")
                        curr_region[n_pt_ind] = reg_idx+1
                        stack.append(n_pt_ind)
                        add_pts_cnt += 1
        output_regions.append(curr_region)
        region_dict = {**curr_region, **region_dict}
        reg_idx +=1
        #print(reg_idx)
    for ind in region_dict:
        #print(ind)
        #inFile.user_data[ind] = region_dict[ind]
        pass
    # for ind in seed_dict:
    #     #print(ind)
    #     inFile.user_data[ind] = 1
    return output_regions

def extract_border_pts(inFile, tree, region_dict):
    border_output = []
    for region in region_dict:
        #print(len(region))
        #print(region)
        if len(region)>10:
            region_coords = []
            ordered_ind_r = []
            for pt_ind in region:
                pt_coord = [inFile.x[pt_ind], inFile.y[pt_ind]]
                region_coords.append(pt_coord)
                ordered_ind_r.append(pt_ind)
            #create a local tree for querying neighbours
            local_tree = spatial.cKDTree(region_coords)
            non_boundary_pts = {}
            for pt_coord in region_coords:
                dist, local_neighbs = local_tree.query(pt_coord, 10)
                #print(local_neighbs)
                local_n_coords=[]
                local_n_ind=[]
                for n_ind in local_neighbs:
                    local_n_coords.append(region_coords[n_ind])
                    local_n_ind.append(ordered_ind_r[n_ind])
                #print(local_n_coords)
                local_border_ind = spatial.ConvexHull(local_n_coords).vertices
                global_border_ind = []
                for local_pt_ind in local_border_ind:
                    global_border_ind.append(local_n_ind[local_pt_ind])
                #print(global_border_ind)
                for pt_ind in local_n_ind:
                    #print(pt_ind)
                    if pt_ind not in global_border_ind:
                        non_boundary_pts[pt_ind]=0
            curr_border =[]
            for pt_ind in region:
                if pt_ind not in non_boundary_pts:
                    curr_border.append(pt_ind)
                    #print(pt_ind)
                    inFile.user_data[pt_ind] = 1
                else:
                    inFile.user_data[pt_ind] = 2
            # for pt_ind in region:
            #     if pt_ind not in non_boundary_pts:
            #         inFile.user_data[ordered_ind_r[pt_ind]] = 1
        else:
            curr_border = []
            for pt_ind in region:
                curr_border.append(pt_ind)
                inFile.user_data[pt_ind] = 1
        border_output.append(curr_border)
    return(border_output)

def border_pts_to_geometry(inFile,border_pt_groups):
    #prepare outfile
    outfile = open("wkt_test.txt", 'w')
    for region_border in border_pt_groups:
        if len(region_border)>4:
            print("starting new region", len(region_border))
            region_border_coord = []
            for border_pt_ind in region_border:
                border_pt_coord = [inFile.x[border_pt_ind], inFile.y[border_pt_ind]]
                region_border_coord.append(border_pt_coord)
            local_tree = spatial.cKDTree(region_border_coord)
            #let's build a graph using knn
            graph_matrix = np.zeros([len(region_border),len(region_border)])
            for pt_coord in region_border_coord:
                dist, local_neighbs = local_tree.query(pt_coord, 5)
                #print(dist,local_neighbs)
                #print(region_border_coord)
                for i in range(1,5):
                    graph_matrix[local_neighbs[0],local_neighbs[i]]=dist[i]
            graph_matrix = np.ma.masked_values(graph_matrix, 0)
            graph = csr_matrix(graph_matrix)
            mst = minimum_spanning_tree(graph).toarray()
            #print(mst)
            out_linestr= "MULTILINESTRING ("
            for orig in range(0,len(region_border)):
                for dest in range(0,len(region_border)):
                    if mst[orig,dest] >0:
                        linestr = '('+str(region_border_coord[orig][0])+' '+ str(region_border_coord[orig][1])+ ',' + str(region_border_coord[dest][0])+ ' '+str(region_border_coord[dest][1])+'),'
                        out_linestr+=linestr
            out_linestr = out_linestr[:-1]+")"
            outfile = open("wkt_test.txt", 'a')
            outfile.write(out_linestr)
            outfile.write("\n")
#code taken from: https://www.sanfoundry.com/python-program-find-minimum-spanning-tree-using-prims-algorithm/

def get_normals_and_merge(roof_surf_list):
    list_normals = []
    idx_to_merge = []
    for roof_surf in roof_surf_list:
        #print(len(roof_surf_list))
        if len(roof_surf)>0:
            normal, ref_pt = fit_plane(roof_surf[0:3])
            list_normals.append(normal)
    for idx in range(0,len(list_normals)-1):
        n_vec = np.array(list_normals[idx])
        for other_idx in range(idx+1,len(list_normals)):
            #print(list_normals)
            other_n_vec = np.array(list_normals[other_idx])
            #dot_product = n_vec[0]*other_n_vec[0] + n_vec[1]*other_n_vec[1] + n_vec[2]*other_n_vec[2]
            dot_product = np.dot(n_vec, other_n_vec)
            if dot_product>0.98:
                intersec_cnt = 0
                for pt in roof_surf_list[idx]:
                    if pt in roof_surf_list[other_idx]:
                        intersec_cnt+=1
                if intersec_cnt>=2:
                    if type(roof_surf_list[idx]) == type(()):
                        roof_surf_list[idx] = [roof_surf_list[idx],roof_surf_list[other_idx]]
                        #print("added")
                        idx_to_merge.append(idx)
                    elif type(roof_surf_list[idx]) == type([]):
                        roof_surf_list[idx].append(roof_surf_list[other_idx])
                        idx_to_merge.append(idx)
                        #print("added2")
                    roof_surf_list[other_idx] = ()
                    #print("erased", other_idx)
                    list_normals[other_idx] = [np.nan, np.nan, np.nan]
            #print(idx, other_idx, dot_product)
    for idx in set(idx_to_merge):
        roof_surf_list[idx] = tuple(merge_geoms(roof_surf_list[idx]))
    roof_surf_list = [x for x in roof_surf_list if len(x)>0]
    #print(list_normals[3][0], type(3.21))
    list_normals = [x for x in list_normals if x[0] is not np.nan]
    #print(list_normals)
    #print(roof_surf_list)
    return roof_surf_list, list_normals

def merge_geoms(geom_list):
    from shapely.ops import cascaded_union
    polygons = []
    for el in geom_list:
        polygons.append(Polygon(el))
    u = cascaded_union(polygons)
    #print(str(u))
    return pgSQL.wkt_to_coords(str(u))

    #print(geom_list)
    #common_pts = 0
    # for geom_idx in range(0,len(geom_list)-1):
    #     for other_geom_idx in range(geom_idx+1,len(geom_list)):
    #         intersec_pts = set(geom_list[geom_idx]).intersection(geom_list[other_geom_idx])
    #         if len(intersec_pts)>=2:
    #             print(geom_idx, other_geom_idx, intersec_pts)
    #             polygon1 = Polygon(geom_list[geom_idx])
    #             polygon2 = Polygon(geom_list[other_geom_idx])
    #             polygons = [polygon1,polygon2]
    #             u = cascaded_union(polygons)
    #             print(u)
                # intersec_idxs = []
                # for inters_pt in intersec_pts:
                #     intersec_idxs.append((geom_list[geom_idx].index(inters_pt), geom_list[other_geom_idx].index(inters_pt)))

            # intersec_list = []
            # for pt_idx in geom_list[geom_idx]:
            #     if pt in geom_list[other_geom_idx]:
            #         intersec_list.append((geom_list[geom_idx].index(pt),geom_list[other_geom_idx].index(pt)))
            #     if len(intersec_list)>2:
            #         print("got sth to merge", intersec_list)


def point_in_polygon(roof_surf_list, ptcloud):
    polygon_list = []
    points_poly_match = []
    for el in roof_surf_list:
        polygon_list.append(Polygon(el))
    for poly_idx in range(0, len(polygon_list)):
        curr_polygon_list = []
        for pt_idx in range(0, len(ptcloud.x)):
            pt = Point(ptcloud.x[pt_idx], ptcloud.y[pt_idx])
            if polygon_list[poly_idx].contains(pt)==True:
                curr_polygon_list.append(pt_idx)
        points_poly_match.append(curr_polygon_list)
    return points_poly_match

t0 = time.time()
inFile, tree = data_h_filter('waalhaven zuid/high_res - Copy.las')
#inFile, tree = data_loader('globalPC/AHN_buildings.las')
deviations = height_deviations(inFile, 0.2)
#print(deviations)
#comp_normals(inFile, tree, 5)
#deviations = comp_height_diff(inFile, tree)
#print(deviations)
#deviations = [0, 1, 2, 3, 4, 6, 7, 13, 14, 16, 18, 19, 20, 21, 28, 30, 31, 32, 33, 35, 36, 43, 44, 45, 46, 48, 49, 50, 55, 56, 58, 59, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 79, 80, 82, 87, 88, 89, 90, 93, 94, 95, 96, 97, 108, 109, 115, 116, 117, 118, 119, 120, 121, 124, 125, 126, 127, 128, 129, 130, 131, 132, 145, 152, 153, 154, 156, 160, 163, 164, 165, 166, 167, 169, 170, 186, 187, 192, 193, 194, 195, 201, 202, 203, 204, 205, 206, 207, 208, 209, 226, 227, 231, 232, 233, 234, 235, 236, 237, 238, 239, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 272, 273, 274, 276, 277, 279, 280, 281, 282, 283, 284, 285, 293, 294, 295, 296, 297, 298, 299, 300, 304, 305, 327, 328, 329, 333, 334, 335, 336, 337, 338, 339, 350, 351, 352, 354, 355, 357, 359, 360, 361, 386, 387, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 414, 415, 416, 417, 418, 419, 421, 423, 424, 426, 427, 455, 456, 461, 463, 464, 466, 484, 485, 488, 489, 496, 497, 498, 527, 528, 529, 533, 534, 535, 536, 537, 539, 559, 560, 603, 607, 608, 609, 610, 611, 612, 636, 637, 638, 641, 646, 647, 648, 683, 684, 689, 690, 691, 692, 693, 694, 720, 721, 722, 723, 724, 725, 726, 727, 728, 765, 766, 767, 769, 770, 771, 773, 774, 775, 776, 805, 806, 807, 808, 809, 810, 849, 850, 852, 855, 857, 859, 860, 889, 892, 893, 894, 896, 897, 899, 900, 901, 908, 909, 910, 923, 924, 925, 943, 944, 946, 947, 948, 954, 988, 989, 990, 991, 992, 993, 994, 995, 996, 1000, 1001, 1002, 1004, 1005, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1042, 1043, 1044, 1081, 1084, 1085, 1086, 1087, 1088, 1091, 1092, 1093, 1101, 1102, 1116, 1117, 1119, 1120, 1122, 1123, 1124, 1144, 1145, 1146, 1147, 1148, 1151, 1153, 1154, 1185, 1188, 1189, 1192, 1194, 1195, 1219, 1220, 1221, 1224, 1226, 1227, 1247, 1248, 1249, 1285, 1288, 1289, 1290, 1315, 1316, 1317, 1320, 1350, 1351, 1352, 1356, 1358, 1391, 1393, 1397, 1399, 1418, 1420, 1421, 1422, 1423, 1425, 1437, 1438, 1456, 1457, 1458, 1459, 1460, 1464, 1496, 1499, 1501, 1503, 1504, 1505, 1527, 1529, 1530, 1531, 1541, 1542, 1543, 1561, 1562, 1563, 1602, 1606, 1607, 1608, 1632, 1633, 1634, 1635, 1643, 1644, 1651, 1669, 1670, 1672, 1675, 1676, 1677, 1710, 1713, 1717, 1718, 1745, 1747, 1782, 1783, 1785, 1789, 1791, 1823, 1824, 1829, 1831, 1893, 1894, 1895, 1901, 1902, 1933, 1934, 1935, 1938, 1939, 1940, 1978, 1979, 1988, 1989, 2006, 2007, 2013, 2052, 2053, 2093, 2094, 2101, 2102, 2103, 2123, 2126, 2127, 2128, 2130, 2164, 2165, 2166, 2167, 2172, 2173, 2209, 2210, 2212, 2214, 2215, 2216, 2220, 2244, 2286, 2287, 2293, 2329, 2330, 2331, 2332, 2333, 2334, 2336, 2340, 2366, 2404, 2405, 2406, 2407, 2413, 2414, 2448, 2449, 2452, 2453, 2454, 2455, 2456, 2458, 2484, 2485, 2488, 2491, 2492, 2525, 2526, 2532, 2533, 2566, 2567, 2569, 2570, 2602, 2603, 2605, 2606, 2607, 2608, 2643, 2644, 2649, 2650, 2719, 2720, 2722, 2724, 2726, 2761, 2762, 2765, 2767, 2768, 2769, 2770, 2838, 2839, 2840, 2879, 2880, 2881, 2882, 2884, 2885, 2886, 2887, 2888, 2889, 2957, 2958, 2959, 2962, 2963, 2980, 2981, 2982, 2983, 2994, 2995, 2998, 3003, 3005, 3006, 3075, 3076, 3093, 3094, 3095, 3096, 3097, 3098, 3113, 3114, 3119, 3120, 3121, 3122, 3189, 3190, 3208, 3209, 3212, 3213, 3214, 3226, 3228, 3231, 3235, 3237, 3238, 3305, 3323, 3324, 3325, 3326, 3329, 3330, 3342, 3344, 3347, 3348, 3349, 3350, 3352, 3354, 3355, 3418, 3419, 3421, 3440, 3441, 3442, 3444, 3445, 3446, 3447, 3448, 3462, 3464, 3465, 3469, 3471, 3472, 3532, 3533, 3534, 3556, 3557, 3560, 3561, 3562, 3563, 3564, 3580, 3583, 3585, 3586, 3587, 3588, 3646, 3647, 3671, 3672, 3673, 3674, 3677, 3678, 3679, 3680, 3694, 3695, 3699, 3700, 3701, 3703, 3704, 3759, 3760, 3761, 3787, 3788, 3789, 3790, 3792, 3793, 3794, 3809, 3813, 3815, 3816, 3868, 3869, 3870, 3899, 3900, 3901, 3903, 3904, 3905, 3906, 3907, 3920, 3921, 3925, 3927, 3928, 3978, 3979, 3980, 4010, 4011, 4012, 4015, 4016, 4017, 4031, 4032, 4033, 4034, 4038, 4040, 4041, 4088, 4089, 4090, 4115, 4116, 4117, 4125, 4126, 4127, 4130, 4131, 4145, 4146, 4147, 4148, 4149, 4153, 4201, 4202, 4239, 4240, 4241, 4242, 4244, 4245, 4246, 4260, 4261, 4262, 4263, 4264, 4265, 4266, 4267, 4268, 4269, 4270, 4274, 4275, 4320, 4321, 4324, 4325, 4361, 4362, 4364, 4365, 4366, 4367, 4368, 4382, 4383, 4384, 4388, 4389, 4390, 4391, 4392, 4394, 4396, 4398, 4399, 4400, 4441, 4442, 4443, 4460, 4461, 4462, 4463, 4470, 4471, 4472, 4484, 4485, 4489, 4490, 4491, 4492, 4505, 4507, 4508, 4509, 4513, 4514, 4515, 4516, 4517, 4518, 4519, 4520, 4521, 4522, 4527, 4566, 4567, 4568, 4580, 4581, 4586, 4589, 4590, 4591, 4592, 4613, 4614, 4615, 4619, 4620, 4621, 4636, 4637, 4646, 4647, 4649, 4651, 4652, 4653, 4654, 4655, 4660, 4697, 4698, 4704, 4705, 4706, 4718, 4719, 4720, 4721, 4723, 4724, 4725, 4746, 4747, 4748, 4749, 4752, 4753, 4754, 4755, 4767, 4768, 4769, 4777, 4778, 4779, 4780, 4781, 4786, 4788, 4791, 4792, 4793, 4827, 4828, 4829, 4839, 4840, 4847, 4848, 4849, 4852, 4853, 4854, 4855, 4875, 4876, 4877, 4878, 4880, 4881, 4882, 4883, 4884, 4898, 4899, 4908, 4909, 4910, 4911, 4912, 4923, 4924, 4925, 4926, 4927, 4957, 4958, 4959, 4960, 4972, 4978, 4979, 4982, 4983, 4984, 4985, 4987, 4988, 4989, 4990, 5009, 5010, 5011, 5012, 5014, 5015, 5016, 5017, 5031, 5032, 5033, 5041, 5042, 5043, 5044, 5045, 5058, 5060, 5089, 5090, 5091, 5095, 5106, 5116, 5117, 5118, 5119, 5122, 5123, 5124, 5134, 5144, 5145, 5146, 5147, 5149, 5150, 5151, 5162, 5163, 5165, 5166, 5167, 5168, 5169, 5170, 5177, 5178, 5179, 5180, 5181, 5182, 5183, 5192, 5194, 5200, 5226, 5227, 5228, 5237, 5238, 5244, 5245, 5254, 5255, 5256, 5257, 5260, 5261, 5262, 5283, 5284, 5285, 5289, 5290, 5303, 5304, 5305, 5306, 5313, 5314, 5316, 5317, 5329, 5331, 5334, 5336, 5337, 5361, 5362, 5363, 5371, 5373, 5381, 5382, 5388, 5389, 5391, 5392, 5393, 5394, 5395, 5397, 5398, 5399, 5400, 5420, 5421, 5422, 5423, 5424, 5425, 5439, 5440, 5444, 5445, 5447, 5462, 5464, 5465, 5466, 5469, 5471, 5472, 5492, 5493, 5494, 5526, 5527, 5528, 5529, 5531, 5532, 5533, 5534, 5535, 5555, 5556, 5570, 5571, 5573, 5574, 5575, 5576, 5578, 5579, 5580, 5600, 5602, 5609, 5610, 5629, 5630, 5658, 5659, 5660, 5663, 5664, 5665, 5668, 5669, 5670, 5707, 5708, 5709, 5710, 5711, 5713, 5715, 5736, 5738, 5741, 5743, 5744, 5745, 5761, 5762, 5763, 5796, 5797, 5798, 5799, 5802, 5803, 5844, 5845, 5846, 5847, 5849, 5871, 5873, 5874, 5875, 5878, 5879, 5880, 5881, 5894, 5895, 5935, 5936, 5937, 5938, 5940, 5941, 5942, 5981, 5982, 5983, 5985, 5987, 5988, 5989, 6014, 6016, 6021, 6022, 6023, 6024, 6034, 6035, 6036, 6056, 6057, 6080, 6081, 6082, 6084, 6085, 6086, 6124, 6125, 6127, 6128, 6157, 6158, 6163, 6164, 6165, 6173, 6174, 6181, 6182, 6201, 6202, 6218, 6219, 6220, 6223, 6224, 6225, 6259, 6260, 6261, 6264, 6265, 6266, 6267, 6268, 6269, 6298, 6299, 6303, 6304, 6305, 6310, 6311, 6312, 6355, 6356, 6357, 6358, 6359, 6360, 6362, 6363, 6364, 6398, 6399, 6400, 6401, 6403, 6436, 6443, 6444, 6448, 6449, 6499, 6500, 6501, 6504, 6505, 6506, 6507, 6534, 6535, 6537, 6538, 6539, 6540, 6576, 6578, 6582, 6583, 6584, 6585, 6586, 6587, 6629, 6637, 6638, 6639, 6641, 6642, 6643, 6644, 6669, 6671, 6672, 6673, 6674, 6712, 6714, 6718, 6719, 6775, 6776, 6777, 6781, 6782, 6783, 6787, 6805, 6806, 6809, 6810, 6811, 6850, 6852, 6855, 6856, 6892, 6893, 6897, 6898, 6913, 6914, 6915, 6916, 6918, 6919, 6920, 6942, 6943, 6945, 6946, 6985, 7024, 7025, 7046, 7047, 7048, 7049, 7050, 7051, 7076, 7077, 7078, 7079, 7080, 7081, 7149, 7150, 7157, 7201, 7202, 7203, 7204, 7205, 7206, 7207, 7260, 7261, 7262, 7284, 7296, 7297, 7328, 7329, 7381, 7382, 7383, 7384, 7404, 7405, 7416, 7427, 7428, 7442, 7443, 7444, 7446, 7495, 7496, 7497, 7499, 7500, 7501, 7502, 7524, 7561, 7562, 7563, 7564, 7565, 7566, 7567, 7613, 7614, 7615, 7616, 7618, 7619, 7620, 7621, 7643, 7678, 7679, 7680, 7681, 7682, 7683, 7684, 7730, 7731, 7732, 7733, 7736, 7737, 7738, 7741, 7742, 7760, 7791, 7792, 7795, 7842, 7843, 7844, 7845, 7848, 7849, 7868, 7900, 7901, 7903, 7950, 7951, 7952, 7953, 7955, 7956, 7957, 7958, 7971, 7972, 8009, 8010, 8011, 8012, 8060, 8061, 8062, 8063, 8066, 8067, 8068, 8117, 8118, 8168, 8169, 8170, 8171, 8173, 8174, 8175, 8176, 8186, 8187, 8222, 8223, 8224, 8273, 8274, 8275, 8276, 8278, 8279, 8280, 8281, 8327, 8328, 8330, 8379, 8380, 8381, 8382, 8384, 8385, 8386, 8387, 8409, 8430, 8431, 8483, 8484, 8485, 8486, 8488, 8489, 8490, 8513, 8532, 8533, 8534, 8535, 8558, 8559, 8560, 8561, 8583, 8584, 8585, 8588, 8589, 8590, 8613, 8631, 8632, 8633, 8634, 8635, 8658, 8659, 8660, 8661, 8662, 8663, 8683, 8684, 8685, 8686, 8688, 8689, 8690, 8691, 8712, 8713, 8729, 8730, 8734, 8735, 8736, 8753, 8754, 8755, 8756, 8757, 8759, 8760, 8782, 8783, 8784, 8785, 8787, 8788, 8789, 8790, 8791, 8811, 8812, 8813, 8828, 8834, 8835, 8836, 8837, 8854, 8855, 8856, 8857, 8858, 8859, 8860, 8861, 8862, 8883, 8884, 8885, 8886, 8888, 8889, 8890, 8891, 8892, 8914, 8927, 8936, 8937, 8956, 8957, 8958, 8959, 8961, 8962, 8963, 8964, 8984, 8985, 8986, 8987, 8990, 8991, 8992, 9001, 9002, 9015, 9024, 9025, 9026, 9037, 9038, 9058, 9059, 9060, 9061, 9064, 9065, 9086, 9087, 9088, 9089, 9091, 9092, 9093, 9094, 9125, 9126, 9137, 9138, 9139, 9158, 9159, 9160, 9163, 9164, 9165, 9186, 9187, 9188, 9189, 9191, 9192, 9193, 9196, 9197, 9199, 9200, 9201, 9216, 9220, 9221, 9222, 9235, 9236, 9255, 9256, 9257, 9258, 9261, 9262, 9263, 9284, 9285, 9286, 9287, 9289, 9290, 9291, 9312, 9313, 9314, 9315, 9332, 9351, 9352, 9353, 9354, 9357, 9358, 9359, 9381, 9382, 9383, 9407, 9408, 9409, 9410, 9428, 9429, 9430, 9450, 9451, 9452, 9453, 9456, 9457, 9479, 9480, 9481, 9482, 9502, 9503, 9507, 9508, 9510, 9531, 9532, 9551, 9552, 9553, 9554, 9556, 9557, 9558, 9605, 9630, 9631, 9651, 9652, 9653, 9654, 9656, 9657, 9658, 9659, 9699, 9700, 9704, 9730, 9749, 9750, 9751, 9752, 9755, 9756, 9757, 9758, 9786, 9787, 9799, 9800, 9801, 9831, 9832, 9833, 9850, 9851, 9852, 9853, 9854, 9857, 9858, 9859, 9898, 9931, 9932, 9933, 9934, 9951, 9952, 9953, 9954, 9956, 9957, 9958, 9959, 9960, 9997, 9998, 9999, 10032, 10033, 10054, 10055, 10058, 10059, 10060, 10061, 10095, 10096, 10097, 10134, 10135, 10156, 10157, 10158, 10162, 10195, 10200, 10201, 10233, 10236, 10237, 10254, 10255, 10256, 10259, 10260, 10261, 10290, 10304, 10305, 10329, 10330, 10331, 10332, 10333, 10349, 10350, 10351, 10352, 10354, 10355, 10356, 10357, 10386, 10403, 10404, 10427, 10428, 10429, 10430, 10449, 10450, 10451, 10452, 10454, 10455, 10483, 10502, 10503, 10504, 10505, 10506, 10507, 10530, 10531, 10551, 10552, 10553, 10554, 10580, 10599, 10600, 10601, 10602, 10603, 10606, 10607, 10629, 10630, 10650, 10651, 10673, 10674, 10675, 10693, 10694, 10695, 10701, 10724, 10725, 10768, 10786, 10787, 10788, 10789, 10793, 10794, 10795, 10796, 10797, 10814, 10815, 10818, 10819, 10820, 10822, 10883, 10884, 10885, 10886, 10892, 10893, 10894, 10920, 10954, 10962, 10976, 10977, 10978, 10979, 10983, 10984, 10985, 11007, 11008, 11012, 11013, 11014, 11015, 11047, 11054, 11073, 11074, 11076, 11077, 11078, 11081, 11082, 11083, 11084, 11103, 11104, 11110, 11167, 11168, 11169, 11173, 11174, 11175, 11176, 11177, 11197, 11205, 11206, 11240, 11246, 11247, 11260, 11261, 11262, 11263, 11264, 11269, 11270, 11271, 11291, 11304, 11305, 11306, 11308, 11309, 11359, 11360, 11361, 11362, 11367, 11368, 11369, 11370, 11391, 11392, 11393, 11408, 11409, 11411, 11457, 11458, 11459, 11460, 11465, 11466, 11467, 11468, 11490, 11491, 11509, 11546, 11547, 11554, 11555, 11560, 11561, 11562, 11563, 11564, 11565, 11608, 11622, 11623, 11632, 11633, 11650, 11651, 11652, 11653, 11658, 11659, 11660, 11705, 11725, 11726, 11742, 11743, 11744, 11745, 11749, 11750, 11751, 11800, 11819, 11820, 11839, 11840, 11841, 11844, 11848, 11849, 11850, 11897, 11898, 11899, 11900, 11928, 11929, 11930, 11935, 11936, 11937, 11938, 11943, 11944, 11945, 11946, 11949, 11950, 11997, 11998, 11999, 12000, 12013, 12034, 12035, 12036, 12037, 12042, 12043, 12044, 12045, 12050, 12098, 12099, 12110, 12130, 12131, 12132, 12133, 12138, 12139, 12140, 12196, 12197, 12198, 12226, 12227, 12228, 12229, 12234, 12235, 12236, 12237, 12295, 12296, 12302, 12303, 12304, 12325, 12326, 12327, 12331, 12332, 12333, 12334, 12335, 12397, 12398, 12401, 12402, 12403, 12430, 12431, 12432, 12433, 12496, 12497, 12498, 12500, 12525, 12526, 12527, 12528, 12551, 12552, 12595, 12596, 12597, 12599, 12622, 12623, 12624, 12625, 12626, 12669, 12670, 12671, 12695, 12696, 12698, 12699, 12700, 12764, 12765, 12766, 12767, 12768, 12769, 12792, 12793, 12794, 12796, 12797, 12830, 12831, 12832, 12833, 12834, 12858, 12859, 12860, 12861, 12865, 12866, 12867, 12888, 12889, 12890, 12891, 12893, 12894, 12926, 12927, 12928, 12929, 12930, 12931, 12933, 12953, 12954, 12955, 12956, 12957, 12961, 12962, 12963, 12983, 12984, 12986, 12987, 12989, 12990, 12991, 12999, 13000, 13028, 13029, 13030, 13031, 13033, 13051, 13052, 13053, 13054, 13058, 13059, 13060, 13061, 13081, 13085, 13086, 13089, 13097, 13098, 13099, 13126, 13127, 13128, 13129, 13130, 13131, 13132, 13148, 13149, 13150, 13151, 13152, 13156, 13157, 13158, 13159, 13179, 13185, 13186, 13196, 13197, 13198, 13199, 13226, 13229, 13244, 13245, 13246, 13247, 13252, 13253, 13254, 13273, 13275, 13276, 13282, 13283, 13284, 13286, 13319, 13320, 13340, 13341, 13342, 13346, 13347, 13348, 13369, 13370, 13379, 13380, 13381, 13382, 13383, 13385, 13436, 13437, 13438, 13439, 13440, 13444, 13445, 13446, 13447, 13467, 13480, 13481, 13483, 13489, 13490, 13535, 13536, 13537, 13538, 13539, 13542, 13543, 13544, 13545, 13546, 13566, 13567, 13581, 13582, 13584, 13585, 13586, 13587, 13588, 13589, 13590, 13591, 13592, 13633, 13634, 13635, 13636, 13637, 13640, 13641, 13642, 13643, 13644, 13660, 13663, 13664, 13665, 13681, 13682, 13683, 13684, 13685, 13686, 13687, 13688, 13690, 13691, 13692, 13733, 13734, 13735, 13736, 13737, 13741, 13742, 13743, 13763, 13783, 13784, 13785, 13786, 13787, 13789, 13790, 13791, 13793, 13795, 13837, 13838, 13839, 13844, 13845, 13846, 13847, 13848, 13867, 13868, 13869, 13889, 13890, 13891, 13892, 13893, 13894, 13895, 13896, 13897, 13898, 13899, 13906, 13907, 13908, 13944, 13945, 13946, 13947, 13948, 13952, 13953, 13954, 13955, 13973, 13974, 13975, 14000, 14001, 14002, 14003, 14004, 14005, 14006, 14007, 14008, 14018, 14054, 14055, 14056, 14057, 14062, 14063, 14064, 14083, 14112, 14113, 14114, 14115, 14117, 14118, 14119, 14120, 14131, 14166, 14167, 14168, 14169, 14173, 14174, 14175, 14176, 14195, 14224, 14225, 14226, 14227, 14228, 14229, 14230, 14231, 14232, 14247, 14248, 14249, 14277, 14278, 14279, 14280, 14281, 14285, 14286, 14287, 14288, 14305, 14306, 14307, 14339, 14340, 14341, 14343, 14344, 14360, 14361, 14362, 14386, 14387, 14388, 14389, 14394, 14395, 14396, 14397, 14416, 14451, 14452, 14453, 14454, 14455, 14456, 14474, 14498, 14499, 14500, 14501, 14506, 14507, 14508, 14509, 14528, 14566, 14567, 14569, 14571, 14591, 14593, 14613, 14614, 14615, 14616, 14621, 14622, 14623, 14642, 14643, 14684, 14685, 14686, 14687, 14704, 14705, 14706, 14727, 14733, 14734, 14735, 14755, 14798, 14800, 14819, 14844, 14847, 14848, 14849, 14850, 14869, 14870, 14889, 14914, 14915, 14916, 14917, 14918, 14934, 14935, 14936, 14959, 14960, 14961, 14962, 14963, 14964, 14965, 14984, 15001, 15002, 15003, 15004, 15026, 15027, 15028, 15029, 15049, 15050, 15074, 15093, 15109, 15110, 15111, 15112, 15137, 15138, 15139, 15140, 15142, 15161, 15205, 15206, 15224, 15225, 15226, 15231, 15254, 15255, 15256, 15257, 15258, 15259, 15277, 15317, 15318, 15319, 15337, 15338, 15339, 15340, 15344, 15345, 15346, 15347, 15348, 15365, 15366, 15367, 15368, 15369, 15370, 15388, 15425, 15426, 15427, 15446, 15447, 15448, 15449, 15454, 15455, 15475, 15478, 15479, 15480, 15483, 15502, 15537, 15555, 15556, 15562, 15563, 15564, 15586, 15587, 15592, 15593, 15594, 15595, 15613, 15614, 15615, 15648, 15666, 15667, 15668, 15669, 15673, 15674, 15675, 15676, 15703, 15704, 15705, 15706, 15707, 15725, 15757, 15758, 15769, 15770, 15776, 15777, 15778, 15783, 15784, 15785, 15815, 15817, 15818, 15820, 15837, 15838, 15866, 15867, 15868, 15886, 15887, 15888, 15889, 15893, 15894, 15895, 15928, 15929, 15948, 15949, 15950, 15977, 15985, 15997, 15998, 15999, 16000, 16005, 16006, 16007, 16008, 16043, 16045, 16058, 16059, 16060, 16062, 16063, 16064, 16086, 16087, 16088, 16107, 16108, 16109, 16114, 16115, 16116, 16117, 16154, 16155, 16158, 16159, 16168, 16169, 16170, 16173, 16174, 16175, 16196, 16216, 16217, 16218, 16219, 16220, 16224, 16225, 16226, 16227, 16265, 16266, 16267, 16268, 16269, 16276, 16278, 16280, 16281, 16300, 16320, 16321, 16322, 16323, 16328, 16329, 16330, 16375, 16376, 16382, 16387, 16403, 16404, 16405, 16425, 16426, 16427, 16432, 16433, 16434, 16480, 16481, 16507, 16508, 16509, 16528, 16529, 16530, 16531, 16535, 16536, 16537, 16584, 16585, 16588, 16589, 16590, 16606, 16625, 16626, 16627, 16628, 16633, 16634, 16635, 16684, 16687, 16688, 16689, 16690, 16699, 16715, 16716, 16719, 16720, 16721, 16722, 16727, 16728, 16729, 16780, 16781, 16783, 16790, 16791, 16792, 16812, 16813, 16814, 16818, 16819, 16820, 16821, 16877, 16879, 16885, 16905, 16906, 16907, 16908, 16912, 16913, 16914, 16915, 16916, 16970, 16971, 16973, 16974, 16976, 16977, 16996, 16997, 16998, 17003, 17004, 17005, 17006, 17062, 17063, 17064, 17066, 17067, 17091, 17094, 17095, 17096, 17156, 17157, 17158, 17159, 17161, 17162, 17185, 17186, 17187, 17188, 17189, 17226, 17227, 17228, 17252, 17254, 17255, 17279, 17280, 17281, 17282, 17321, 17322, 17323, 17324, 17347, 17348, 17350, 17414, 17415, 17416, 17417, 17418, 17420, 17421, 17444, 17445, 17447, 17448, 17508, 17509, 17510, 17511, 17515, 17516, 17539, 17540, 17541, 17542, 17544, 17601, 17602, 17603, 17604, 17610, 17611, 17634, 17635, 17636, 17637, 17638, 17640, 17696, 17697, 17698, 17699, 17700, 17704, 17705, 17706, 17707, 17728, 17731, 17732, 17734, 17736, 17737, 17787, 17788, 17789, 17790, 17795, 17796, 17797, 17798, 17819, 17825, 17826, 17828, 17830, 17831, 17879, 17880, 17881, 17882, 17883, 17888, 17889, 17890, 17909, 17911, 17920, 17922, 17972, 17973, 17974, 17978, 17979, 17980, 17981, 17982, 18001, 18002, 18003, 18015, 18017, 18018, 18066, 18067, 18068, 18073, 18074, 18075, 18076, 18096, 18097, 18108, 18109, 18110, 18111, 18156, 18157, 18158, 18159, 18160, 18164, 18165, 18166, 18167, 18188, 18203, 18204, 18205, 18206, 18207, 18208, 18209, 18233, 18252, 18253, 18254, 18255, 18260, 18261, 18262, 18263, 18283, 18284, 18285, 18299, 18300, 18301, 18302, 18303, 18304, 18305, 18307, 18308, 18309, 18310, 18311, 18335, 18336, 18352, 18353, 18354, 18355, 18360, 18361, 18362, 18381, 18382, 18383, 18401, 18402, 18404, 18405, 18406, 18407, 18408, 18410, 18413, 18414, 18451, 18452, 18453, 18454, 18460, 18461, 18462, 18482, 18503, 18504, 18505, 18506, 18507, 18508, 18509, 18510, 18511, 18550, 18551, 18552, 18553, 18558, 18559, 18560, 18579, 18580, 18581, 18587, 18588, 18604, 18605, 18607, 18608, 18610, 18611, 18612, 18621, 18622, 18651, 18652, 18653, 18654, 18655, 18659, 18660, 18661, 18681, 18682, 18710, 18711, 18713, 18714, 18715, 18716, 18717, 18719, 18731, 18757, 18758, 18759, 18760, 18761, 18765, 18766, 18767, 18768, 18788, 18814, 18815, 18816, 18817, 18818, 18819, 18820, 18821, 18822, 18836, 18861, 18862, 18863, 18864, 18869, 18870, 18871, 18872, 18893, 18894, 18895, 18921, 18923, 18924, 18925, 18948, 18949, 18968, 18969, 18970, 18971, 18975, 18976, 18977, 18978, 18998, 18999, 19000, 19033, 19034, 19035, 19038, 19055, 19056, 19057, 19059, 19060, 19076, 19077, 19078, 19079, 19084, 19085, 19086, 19106, 19107, 19108, 19145, 19146, 19148, 19149, 19150, 19167, 19168, 19169, 19190, 19191, 19192, 19196, 19197, 19198, 19199, 19219, 19257, 19258, 19261, 19278, 19279, 19280, 19281, 19305, 19310, 19311, 19330, 19331, 19332, 19370, 19371, 19390, 19416, 19417, 19418, 19419, 19440, 19441, 19453, 19482, 19499, 19500, 19525, 19526, 19548, 19549, 19561, 19562, 19563, 19564, 19565, 19589, 19590, 19591, 19592, 19593, 19594, 19613, 19655, 19656, 19657, 19671, 19672, 19673, 19698, 19699, 19700, 19701, 19702, 19703, 19704, 19723, 19767, 19768, 19781, 19782, 19783, 19784, 19785, 19789, 19809, 19810, 19811, 19812, 19831, 19872, 19887, 19888, 19889, 19893, 19894, 19895, 19896, 19916, 19917, 19919, 19920, 19921, 19939, 19940, 19978, 19993, 19994, 19995, 19996, 20001, 20002, 20003, 20023, 20024, 20044, 20045, 20081, 20082, 20083, 20096, 20097, 20099, 20100, 20103, 20104, 20105, 20132, 20133, 20134, 20153, 20187, 20188, 20200, 20201, 20202, 20208, 20209, 20210, 20211, 20237, 20238, 20239, 20241, 20242, 20243, 20260, 20261, 20293, 20294, 20308, 20309, 20310, 20316, 20317, 20318, 20345, 20346, 20347, 20348, 20350, 20352, 20368, 20370, 20405, 20406, 20413, 20414, 20415, 20416, 20417, 20423, 20424, 20445, 20457, 20458, 20459, 20460, 20461, 20476, 20477, 20478, 20479, 20506, 20507, 20520, 20521, 20522, 20527, 20528, 20529, 20564, 20565, 20566, 20581, 20582, 20583, 20585, 20610, 20611, 20623, 20624, 20625, 20626, 20627, 20632, 20633, 20670, 20672, 20673, 20674, 20684, 20685, 20686, 20687, 20688, 20690, 20711, 20712, 20713, 20726, 20727, 20728, 20729, 20734, 20735, 20736, 20737, 20775, 20776, 20777, 20778, 20779, 20780, 20781, 20788, 20789, 20790, 20791, 20793, 20794, 20812, 20813, 20814, 20815, 20816, 20828, 20829, 20836, 20837, 20838, 20839, 20882, 20883, 20884, 20885, 20893, 20894, 20895, 20897, 20914, 20915, 20924, 20925, 20926, 20929, 20930, 20931, 20936, 20937, 20938, 20983, 20984, 20985, 20987, 20994, 20995, 20997, 21012, 21023, 21027, 21028, 21029, 21030, 21035, 21036, 21088, 21089, 21092, 21093, 21094, 21097, 21098, 21109, 21110, 21111, 21126, 21127, 21128, 21129, 21134, 21135, 21136, 21137, 21187, 21188, 21189, 21190, 21192, 21194, 21195, 21205, 21221, 21222, 21223, 21229, 21230, 21278, 21279, 21280, 21281, 21282, 21284, 21292, 21309, 21310, 21311, 21317, 21318, 21374, 21378, 21382, 21383, 21384, 21394, 21395, 21401, 21402, 21403, 21404, 21405, 21409, 21410, 21411, 21412, 21443, 21466, 21467, 21469, 21470, 21471, 21473, 21474, 21492, 21493, 21494, 21495, 21498, 21499, 21500, 21501, 21502, 21521, 21560, 21561, 21563, 21564, 21565, 21566, 21588, 21591, 21592, 21593, 21594, 21595, 21654, 21655, 21656, 21657, 21658, 21659, 21660, 21684, 21685, 21686, 21687, 21726, 21727, 21728, 21752, 21753, 21755, 21778, 21779, 21819, 21820, 21821, 21822, 21845, 21846, 21848, 21849, 21907, 21908, 21909, 21913, 21936, 21937, 21938, 21939, 21941, 22000, 22001, 22002, 22007, 22008, 22009, 22032, 22035, 22036, 22092, 22093, 22094, 22095, 22096, 22101, 22102, 22122, 22123, 22124, 22126, 22127, 22187, 22188, 22189, 22190, 22191, 22194, 22195, 22196, 22197, 22218, 22219, 22223, 22224, 22225, 22227, 22262, 22263, 22283, 22284, 22285, 22291, 22292, 22293, 22294, 22313, 22314, 22315, 22320, 22321, 22322, 22323, 22324, 22352, 22353, 22375, 22376, 22377, 22378, 22379, 22383, 22384, 22385, 22386, 22405, 22406, 22407, 22408, 22415, 22416, 22418, 22419, 22420, 22421, 22471, 22472, 22473, 22479, 22480, 22481, 22501, 22514, 22515, 22516, 22517, 22518, 22564, 22565, 22566, 22567, 22568, 22572, 22573, 22574, 22575, 22594, 22595, 22610, 22611, 22612, 22613, 22614, 22615, 22643, 22644, 22659, 22660, 22661, 22662, 22666, 22667, 22668, 22669, 22689, 22707, 22709, 22710, 22711, 22752, 22753, 22754, 22755, 22756, 22761, 22762, 22763, 22783, 22800, 22801, 22802, 22803, 22804, 22806, 22807, 22808, 22809, 22810, 22846, 22847, 22848, 22849, 22850, 22854, 22855, 22856, 22857, 22879, 22899, 22900, 22901, 22902, 22903, 22906, 22908, 22910, 22945, 22946, 22947, 22948, 22953, 22954, 22955, 22975, 22976, 22997, 22998, 22999, 23000, 23001, 23002, 23009, 23011, 23042, 23043, 23044, 23045, 23046, 23050, 23051, 23052, 23053, 23073, 23099, 23100, 23101, 23102, 23103, 23113, 23114, 23140, 23141, 23142, 23143, 23147, 23148, 23149, 23171, 23180, 23181, 23182, 23200, 23201, 23202, 23203, 23204, 23207, 23219, 23220, 23221, 23245, 23246, 23247, 23248, 23254, 23255, 23274, 23275, 23276, 23306, 23307, 23308, 23309, 23310, 23311, 23312, 23313, 23328, 23333, 23349, 23350, 23351, 23352, 23357, 23358, 23359, 23379, 23380, 23409, 23410, 23411, 23412, 23413, 23415, 23432, 23436, 23452, 23453, 23454, 23455, 23456, 23460, 23461, 23462, 23481, 23482, 23483, 23516, 23517, 23518, 23519, 23520, 23521, 23523, 23525, 23543, 23545, 23561, 23562, 23563, 23564, 23569, 23570, 23571, 23591, 23629, 23630, 23632, 23651, 23652, 23653, 23654, 23669, 23670, 23671, 23672, 23674, 23675, 23676, 23677, 23678, 23679, 23699, 23700, 23701, 23742, 23745, 23746, 23747, 23763, 23764, 23765, 23766, 23788, 23789, 23790, 23811, 23854, 23872, 23873, 23875, 23895, 23896, 23897, 23898, 23899, 23900, 23919, 23920, 23921, 23936, 23937, 23938, 23939, 23962, 23963, 23964, 23965, 23966, 23969, 23978, 23979, 23980, 23986, 23987, 23988, 24010, 24011, 24031, 24040, 24041, 24042, 24043, 24068, 24069, 24070, 24083, 24084, 24085, 24090, 24091, 24135, 24143, 24144, 24145, 24149, 24171, 24172, 24173, 24174, 24175, 24176, 24186, 24187, 24195, 24196, 24197, 24238, 24239, 24246, 24247, 24248, 24249, 24252, 24253, 24254, 24255, 24276, 24277, 24278, 24279, 24280, 24299, 24340, 24347, 24348, 24349, 24354, 24355, 24356, 24380, 24381, 24382, 24383, 24384, 24385, 24404, 24442, 24443, 24444, 24450, 24451, 24452, 24457, 24458, 24459, 24460, 24489, 24491, 24511, 24545, 24546, 24547, 24556, 24557, 24558, 24562, 24563, 24564, 24565, 24598, 24599, 24601, 24620, 24621, 24654, 24655, 24656, 24667, 24668, 24669, 24673, 24674, 24675, 24707, 24708, 24710, 24711, 24729, 24730, 24761, 24771, 24772, 24773, 24774, 24778, 24779, 24780, 24781, 24802, 24803, 24813, 24814, 24815, 24817, 24819, 24835, 24836, 24837, 24865, 24878, 24879, 24884, 24885, 24886, 24907, 24908, 24922, 24923, 24926, 24939, 24940, 24941, 24943, 24944, 24968, 24969, 24970, 24982, 24983, 24984, 24988, 24989, 24990, 24991, 25030, 25031, 25032, 25033, 25034, 25048, 25049, 25050, 25051, 25074, 25087, 25088, 25089, 25090, 25093, 25094, 25095, 25137, 25138, 25139, 25140, 25152, 25154, 25155, 25177, 25190, 25191, 25196, 25197, 25198, 25199, 25240, 25241, 25242, 25243, 25244, 25245, 25246, 25247, 25254, 25255, 25256, 25258, 25259, 25260, 25279, 25280, 25281, 25290, 25291, 25293, 25294, 25295, 25296, 25298, 25299, 25300, 25301, 25345, 25346, 25347, 25348, 25349, 25350, 25356, 25357, 25359, 25360, 25375, 25376, 25377, 25391, 25392, 25393, 25397, 25398, 25399, 25400, 25401, 25446, 25448, 25450, 25454, 25455, 25456, 25457, 25458, 25472, 25489, 25490, 25491, 25492, 25496, 25497, 25546, 25547, 25548, 25550, 25551, 25552, 25553, 25555, 25556, 25557, 25568, 25585, 25586, 25591, 25592, 25593, 25594, 25595, 25644, 25645, 25647, 25651, 25659, 25676, 25677, 25678, 25679, 25680, 25681, 25685, 25686, 25687, 25741, 25742, 25743, 25744, 25750, 25751, 25758, 25759, 25768, 25769, 25770, 25774, 25775, 25776, 25777, 25832, 25833, 25834, 25835, 25848, 25849, 25859, 25860, 25861, 25862, 25863, 25865, 25866, 25867, 25868, 25870, 25926, 25927, 25928, 25929, 25930, 25931, 25932, 25953, 25954, 25959, 25960, 25961, 25962, 26021, 26022, 26024, 26025, 26047, 26048, 26049, 26050, 26051, 26052, 26053, 26054, 26092, 26093, 26116, 26117, 26121, 26144, 26145, 26185, 26186, 26187, 26210, 26211, 26213, 26214, 26215, 26276, 26277, 26278, 26279, 26304, 26305, 26306, 26307, 26308, 26309, 26310, 26370, 26371, 26372, 26373, 26377, 26400, 26401, 26402, 26403, 26404, 26405, 26463, 26464, 26465, 26466, 26467, 26470, 26471, 26472, 26493, 26494, 26495, 26496, 26554, 26555, 26556, 26561, 26562, 26563, 26564, 26583, 26584, 26585, 26586, 26587, 26644, 26645, 26646, 26651, 26652, 26653, 26674, 26675, 26676, 26677, 26728, 26729, 26730, 26731, 26732, 26736, 26737, 26738, 26739, 26740, 26741, 26759, 26760, 26761, 26762, 26764, 26765, 26817, 26818, 26819, 26820, 26824, 26825, 26826, 26827, 26845, 26846, 26847, 26849, 26877, 26878, 26879, 26880, 26898, 26899, 26900, 26901, 26905, 26906, 26907, 26926, 26927, 26928, 26929, 26930, 26959, 26960, 26977, 26978, 26979, 26980, 26981, 26984, 26985, 26986, 26987, 27006, 27007, 27009, 27010, 27011, 27012, 27013, 27014, 27058, 27059, 27060, 27065, 27066, 27067, 27068, 27087, 27088, 27089, 27090, 27091, 27092, 27093, 27094, 27095, 27096, 27097, 27098, 27136, 27137, 27138, 27139, 27140, 27143, 27144, 27145, 27146, 27164, 27165, 27166, 27167, 27168, 27169, 27170, 27171, 27172, 27173, 27174, 27178, 27212, 27213, 27214, 27215, 27221, 27222, 27240, 27241, 27242, 27243, 27244, 27245, 27246, 27247, 27255, 27287, 27288, 27289, 27290, 27296, 27297, 27313, 27314, 27315, 27316, 27317, 27318, 27319, 27320, 27323, 27330, 27340, 27357, 27358, 27359, 27360, 27361, 27364, 27365, 27366, 27367, 27368, 27387, 27388, 27389, 27390, 27391, 27392, 27393, 27394, 27395, 27404, 27405, 27406, 27407, 27430, 27431, 27432, 27433, 27436, 27437, 27438, 27439, 27440, 27441, 27456, 27457, 27458, 27459, 27460, 27461, 27462, 27463, 27464, 27477, 27478, 27498, 27499, 27500, 27501, 27502, 27505, 27506, 27507, 27523, 27524, 27525, 27526, 27527, 27528, 27529, 27530, 27531, 27546, 27568, 27569, 27570, 27571, 27576, 27577, 27578, 27594, 27595, 27596, 27597, 27617, 27618, 27619, 27635, 27636, 27637, 27638, 27642, 27643, 27644, 27661, 27662, 27663, 27664, 27665, 27667, 27685, 27686, 27702, 27703, 27704, 27705, 27706, 27711, 27712, 27713, 27729, 27730, 27731, 27732, 27752, 27753, 27754, 27771, 27775, 27776, 27777, 27778, 27779, 27796, 27797, 27798, 27799, 27800, 27801, 27820, 27821, 27822, 27844, 27845, 27846, 27847, 27866, 27867, 27869, 27870, 27890, 27891, 27910, 27911, 27912, 27913, 27914, 27933, 27934, 27935, 27936, 27938, 27939, 27940, 27958, 27959, 28002, 28003, 28004, 28005, 28008, 28010, 28028, 28071, 28072, 28073, 28134, 28135, 28136, 28137, 28138, 28139, 28140, 28159, 28160, 28199, 28200, 28201, 28202, 28203, 28204, 28222, 28258, 28259, 28260, 28261, 28263, 28265, 28281, 28282, 28283, 28316, 28317, 28318, 28319, 28320, 28322, 28339, 28372, 28373, 28374, 28375, 28376, 28377, 28378, 28380, 28397, 28398, 28399, 28428, 28429, 28430, 28431, 28432, 28433, 28434, 28435, 28436, 28450, 28451, 28452, 28453, 28454, 28455, 28480, 28481, 28482, 28483, 28484, 28485, 28486, 28499, 28501, 28503, 28504, 28505, 28506, 28528, 28529, 28530, 28531, 28532, 28533, 28545, 28546, 28550, 28551, 28552, 28553, 28575, 28576, 28577, 28578, 28579, 28580, 28581, 28590, 28591, 28592, 28597, 28598, 28599, 28619, 28620, 28621, 28622, 28623, 28624, 28626, 28633, 28635, 28638, 28639, 28640, 28657, 28658, 28659, 28660, 28662, 28667, 28668, 28669, 28672, 28674, 28675, 28676, 28691, 28692, 28693, 28695, 28697, 28698, 28699, 28701, 28702, 28704, 28705, 28718, 28719, 28721, 28722, 28728, 28729, 28731, 28732, 28733, 28742, 28743, 28744, 28746, 28747, 28748, 28751, 28757, 28758, 28759, 28764, 28768, 28770, 28771, 28772, 28773, 28774, 28778, 28780, 28781, 28782, 28783, 28784]
#regions = grow_regions(deviations, tree, inFile, 5)
#regions = [{8435: 28, 8433: 28, 8332: 28, 8334: 28, 8434: 28, 8436: 28, 8330: 28, 8333: 28, 8227: 28, 8122: 28, 8014: 28, 8226: 28, 8225: 28, 8331: 28, 8121: 28, 8120: 28, 8013: 28, 7904: 28, 8012: 28, 7903: 28}]
#border_pt_groups = extract_border_pts(inFile,tree, regions)
#print(border_pt_groups)
#border_pt_groups = [[65, 198, 286, 401, 468, 541, 614, 778, 861, 1045, 948, 393, 336, 233, 152, 87, 535, 689, 92, 162, 244, 415, 483, 635, 719, 984, 1184, 987, 1083, 1496, 1498, 1934, 2162, 2283, 2403, 2405, 2995, 2997, 3110, 3227, 3230, 3692, 3463, 3579, 4259, 4261, 4382, 4504, 4633, 5030, 5165, 5303, 5304], [2726, 2608, 2724], [3563, 3796, 4247, 4755, 4884, 5009, 4875, 4484, 4239, 3899, 3671, 3323, 3093, 3098, 2983, 2981, 2979, 5420, 5555, 5556, 5425, 5290, 5152], [4653, 4914, 5183, 5179, 5448, 5446, 5583, 5717, 5990, 6270, 6405, 6542, 6676, 7080, 7331, 7446, 7682, 7684, 7795, 5711, 5983, 6125, 5577, 5444, 5167, 5031, 4898, 4636, 5034, 4770, 4383, 4262, 4031, 4645, 4647, 5177, 4268, 3920, 6399, 6672, 6809, 7078, 7204, 4396], [4718, 4847, 4585, 4459, 4338, 4339, 4463, 4726, 5116, 5391, 5125, 5796, 6079, 5671, 5804, 6499, 6226, 6365, 6913, 6775, 7048, 7179, 7180, 7181, 7054, 6783], [7286, 7522, 7641, 7758, 7866, 7406, 7645, 8088, 7980, 8200, 8302, 8512, 8612, 8411, 8306, 8514, 8613, 7157], [7384, 7621, 7261, 7260, 7258, 7380, 7613, 7730, 8060, 7738, 7958, 8176, 8378, 8482, 8281, 8387, 8782, 8984, 8692, 8891, 9094, 9186, 9479, 9483, 9579, 9385, 9291, 8590], [8435, 8433, 8334, 8436, 8330, 8227, 7904, 8012, 7903], [8813, 8812, 9013, 9115, 9017, 9509, 9602, 9802, 9895, 10098, 10193, 10387, 10579, 10770, 11726, 9411, 9122, 9316, 9127, 8928, 8825, 8630, 8730, 8636, 9137, 9139, 9236, 9331, 9428, 9530, 9631, 9830, 10232, 10235, 10531, 10528, 10628, 10818, 10820, 11012, 11304, 11112, 11407, 11409, 11704, 11899, 11996, 12097, 12396, 12297, 12496, 12794, 12987, 13087, 13284, 13379, 13581, 13583, 13681, 13892, 14110, 14226, 14338, 14450, 14565, 14683, 14912, 15136, 15480, 15478, 15927, 15929, 16479, 16481, 16877, 16971, 17158, 17251, 17821, 18000, 18191, 18185, 18384, 18483, 18679, 19110, 19106, 19768, 19654, 19976, 20714, 20608, 20815, 20710, 21111, 21203, 21756, 22129, 22419, 22519, 22712, 22805, 22811, 22912, 23766, 23638, 23860, 24497, 24285, 24716, 25358, 25257, 25554, 25647, 25549, 25449, 25351, 25141, 24926, 24816, 24709, 24603, 24386, 24600, 24490, 23969, 23853, 23747, 23744, 23631, 23416, 23414, 23522, 23003, 23529, 25145, 25038, 24838, 24945, 25459, 25457, 25649, 25743, 25651, 25835, 25836, 24622, 24301, 23545, 23654, 23333, 23204, 23412, 23201, 23964, 24923, 25139, 24174, 25548, 22128, 21376, 21096, 20792, 20372, 20046, 19725, 19392, 19169, 19057, 19044, 18823, 18509, 18415, 18312, 18838, 18306, 17925, 17737, 17734, 17828, 17542, 17545, 17351, 17159, 16594, 16493, 16490, 16173, 16064, 15951, 15615, 15502, 15389, 15162, 15051, 14476, 14249, 14132, 13908, 13692, 15268, 15040, 14810, 14695, 14465, 14117, 14342, 14568, 14915, 14918, 15032, 15481, 15483, 15816, 15820, 15932, 16046, 16157, 16159, 16482, 16586, 16685, 16686, 16687, 16591, 16385, 16171, 15828, 15491, 13788, 13487, 13383, 13190, 13187, 12891, 12600, 12400, 12399, 12201, 12298, 11999, 11902, 11706, 11412, 11410, 11210, 11115, 11113, 10824, 10535, 10821, 10532, 10432, 10429, 10138, 10036, 10034, 9731, 9733, 9632, 9239, 9140, 9041, 8836, 8736, 8735, 14708, 16887, 16697, 16403, 16508, 16194, 15870, 15645, 15425, 15091, 15321, 15095, 14985, 14758, 14526, 14308, 14081, 13765, 13564, 13465, 13372, 13177, 13079, 13660, 14752, 18927, 18716, 18606, 17922, 19034, 18924, 19591, 19917, 20027, 19813, 19918, 19705, 19482, 19479, 19262, 19036, 20240, 20352, 20567, 20671, 20781, 20986, 20989, 21192, 21095, 20895, 20680, 20791, 20687, 20249, 20033, 19927, 19601, 19379, 19156, 20358, 20573, 19259, 20884, 20777, 21189, 20348, 20691, 20795, 20998, 20995, 21285, 21939, 21657, 18895, 17637, 17733, 18013, 17921, 18204, 18400, 18403, 18921, 19031, 18816, 19143, 19475, 19368, 19808, 19916, 19810, 19700, 15139, 8327, 8222, 8119, 8009, 7899, 7902, 7560, 7681, 7325, 7445, 7074, 6944, 6805, 6808, 6533, 6536, 6395, 6258, 6262, 5982, 5980], [9157, 9351, 9458, 9659, 9758, 9951, 9850, 10253, 10349, 10650, 10651, 10357, 10456, 10061, 8964, 8862, 8662, 9057, 8854, 8656, 8556, 8460, 8461], [10404, 10406, 10607, 10797, 10987, 10692, 10882, 11468, 11260, 12045, 12237, 11839, 11650, 11554, 11935, 12530, 12626, 12326, 12130, 12226, 12622, 11853, 12052, 12149, 12148, 11272], [11492, 11493, 11591, 11592, 11490], [12963, 13061, 12866, 13159, 13255, 13050, 13447, 13644, 13338, 13435, 13848, 13955, 14176, 14288, 14509, 14385, 14497, 14612, 14842, 14964, 15074, 15073, 14736, 14053, 13732, 13944, 13147, 12764, 12673, 12571], [13489, 13490, 13591, 13590, 13592], [14889, 14891, 15109, 15677, 15446, 15666, 15897, 15997, 16320, 16636, 16823, 16719, 17098, 17282, 17090, 17184, 16997, 16905, 17278, 17375, 16424, 16436, 16227, 15457, 15347, 15232], [17414, 17707, 17612, 17517, 17230, 17228, 17227, 17226, 17982, 17601, 17695, 17879, 18262, 18167, 18462, 18768, 18872, 18549, 18757, 19086, 19311, 19419, 18967, 19076, 19190, 19413, 19523, 19526, 18451, 18065, 18156, 18251], [19455, 19453, 19561, 19781, 20200, 20415, 20623, 20634, 20426, 20107, 20003, 19790, 20839, 21038, 20726, 20828, 21126, 21503, 21595, 21871, 21682, 21777, 21493, 21309, 21320, 21232], [21656, 21465, 21375, 21278, 21186, 21087, 20881, 20774, 20669, 20562, 20345, 20347, 20239, 20238, 22034, 22415, 22417, 22800, 22802, 22997, 22998, 22408, 22215, 22311, 22596, 22499, 22691, 22687, 22977, 22876, 22973, 23074, 23274, 23481, 23699, 23919, 24136, 24240, 24133, 24548, 24545, 24762, 24759, 24968, 25470, 25569, 25752, 25656, 25748, 25840, 25841], [22225], [22563, 22857, 23053, 23359, 23571, 23348, 23452, 23561, 23670, 23784, 24011, 24010, 23895, 23790, 23679, 23139, 22944, 22658, 22282, 22186, 22481, 22386, 22198, 22009, 21727, 21725, 21906, 21999, 21728, 21729, 21914, 22669], [24566, 24347, 24556, 24771, 24982, 25087, 25401, 25199, 24992, 24781, 25595, 25778, 25585, 25391, 25677, 25953, 26146, 26147, 25962, 26143, 26047, 24246, 23936, 23829, 23828, 24045, 24357], [26472, 26283, 26184, 26092, 26826, 26987, 26897, 27057, 27223, 27146, 27368, 27441, 27285, 27712, 27634, 27702, 27770, 27909, 27982, 27983, 27847, 27497, 27578, 27135, 26816, 26462, 26369, 26642, 26739, 26653, 26564, 26093, 26188], [26763, 26766, 26850, 27548, 27332, 27257, 27821, 27893, 27823, 27753, 28161, 28284, 28670, 28726, 28723, 28661, 28593, 28501, 28400, 28342, 28262, 28007, 27937, 27598, 27666, 26848, 27091, 27168, 27526, 27527, 28375, 28431, 28721, 28693, 28579, 28532, 28377, 27664, 27597, 27462, 26496, 26118, 26023, 25928, 25259, 25154, 25049, 25258, 26499, 26121, 26402], [638], [4960, 5093, 4962, 5232, 5097, 5235, 5504, 5503], [5335, 5468, 5607, 5742, 5877, 6020, 6582], [7564], [12735, 12732, 12925, 12928, 12930, 12931, 12835], [12903, 12904, 12997, 13001, 13197, 13198, 13195, 13199], [13031, 13032, 13130, 13132, 13229, 13129], [13485, 13188, 13483, 14006, 14114, 14115, 14453, 14341, 14005, 15140, 14684, 15818, 16156, 13684], [13889, 13890], [21282], [24069, 24171, 24066, 23962, 23963, 24382, 24485, 24488, 24597, 24705, 24812, 24920, 25028, 25138, 25242, 25447, 25444, 25643, 25831, 25833, 26020, 26115, 26209, 26306, 26495, 26583, 26761, 26847, 26758, 26844, 27008, 27005, 27243, 27313, 27458, 27386, 27525, 27593, 27663, 27798, 27865, 28002, 28200, 28197, 28372, 28430, 28528, 28621, 28657, 28720, 28717, 28757, 28779, 28777, 28783], [25929, 26024], [28550, 28596, 28505, 28554, 28504, 28674, 28731, 28706, 28749, 28762, 28751, 28774, 28780]]
#border_pts_to_geometry(inFile,border_pt_groups)
map_curvature(inFile, tree, 20)

# inFile, tree = data_loader('globalPC/AHN_buildings.las')
# test = pgSQL.CityDB_connection("testdb", "127.0.0.1")
# bbox = test.get_bbox(232)
# reduced_pc, reduced_tree = crop_pc(inFile,bbox,232)
# roof_surf_list = test.get_roof_geom(232)
# roof_surf, n_vecs = get_normals_and_merge(roof_surf_list)
# point_in_poly_list = point_in_polygon(roof_surf, reduced_pc)
# dev_seeds = height_deviations_vertical(reduced_pc, 0.4, point_in_poly_list, roof_surf, n_vecs)
# for idx in roof_surf:
#     if dev_seeds[idx] != []:
#         local_coord_list = []
#         for pt in point_in_poly_list[idx]:
#             local_coord_list.append([reduced_pc.x[pt],reduced_pc.y[pt],reduced_pc.z[pt]])
#         local_tree = spatial.cKDTree(local_coord_list)
#         regions = grow_regions(dev_seeds[idx], local_tree, reduced_pc, 5)
#         border_pt_groups = extract_border_pts(inFile, local_tree, regions)

#loc_idx = 0
# for id in ids:
#     reduced_pc, reduced_tree = crop_pc(inFile,bboxs[loc_idx],id)
#     roof_surf = test.get_roof_geom(id)
#     loc_idx +=1
#     print(roof_surf)
t1 = time.time()
print("comp time", t1-t0)















# inFile = File('kerk/kerk_bldg - copy.las', mode= 'rw')
# tree = spatial.KDTree(list(zip(inFile.X, inFile.Y)))
# for index in range(0,len(inFile.X)):
#     neighbors = (tree.query_ball_point([inFile.X[index], inFile.Y[index]], 500))
#     valuelist = []
#     for pt in neighbors:
#         valuelist.append(inFile.intensity[pt]*100.0)
#     inFile.Z[index] = sum(valuelist)/len(valuelist)
# def extend_patch(index, radius=500, ref_intensity, tolerance):
#     neighbors = (tree.query_ball_point([inFile.X[index], inFile.Y[index]],radius))
#     patch_list = []
#     for pt in neighbors:
#         intensity_n = inFile.Intensity[pt])
#         if abs(intensity_n-ref_intensity)<tolerance:

# def covMat(points):
#     sumx=0
#     sumy=0
#     sumz=0
#     covmat=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
#     for pt in points:
#         sumx=sumx+pt.X
#         sumy=sumy+pt.Y
#         sumz=sumz+pt.Z
#     cx=sumx/len(points)
#     cy=sumy/len(points)
#     cz=sumz/len(points)
#     covLst=[]
#     for pt in points:
#         ppbar=np.matrix([[pt.X-cx,pt.Y-cy,pt.Z-cz]])
#         ppbart=ppbar.transpose()
#         Mat3=np.multiply(ppbart,ppbar)
#         covmat=covmat+Mat3
#     c=(1/float(len(points)))*covmat
#     return c
#
# Nvector=[]
# sigma=[]
# for item in points:
#     a = np.linalg.eig(covMat(item))
#     eigenvalues=list((a[0]))
#     eigenvectors=(a[1]).tolist()
#     lamda0=min(eigenvalues)
#     index0=eigenvalues.index(lamda0)
#     n=(rg.Vector3d(eigenvectors[index0][0],eigenvectors[index0][1],eigenvectors[index0][2]))
#     sigma.append(lamda0/float(sum(eigenvalues)))
#     Nvector.append(n)
