# import plotly.plotly as py
# import plotly.graph_objs as go
# import numpy as np

#test = go.mesh3d(x=[0,0,1,1,2,2,3,3], y=[0,2,2,1,1,2,2,0], z=[0,0,0,0,0,0,0,0], alphahull=5, opacity=0.4, color='#00FFFF')
#see: https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html
#code taken from: https://gist.github.com/dwyerk/10561690

from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Delaunay, Voronoi
import numpy as np
import math
import shapely.geometry as geometry
import shapely.wkt as wkt

def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set of points.

    @param points: Iterable container of points.
    @param alpha: alpha value to influence the gooeyness of the border. Smaller
                  numbers don't fall inward as much as larger numbers. Too large,
                  and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense in computing an alpha
        # shape.
        return geometry.MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])

    coords = np.array([point.coords[0] for point in points])

    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    print(tri.vertices)
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        # Semiperimeter of triangle
        s = (a + b + c) / 2.0

        # Area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)

        # Here's the radius filter.
        # print circum_r
        if circum_r**2 < alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
        else:
            print("refused",pa,pb,pc, circum_r)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

def voronoi_shape(all_points, region_dict, distance_list, bounding_polygon):
    #print(geometry.MultiPoint(all_points))
    #print("input",all_points,",", region_dict,",", bounding_geometry,",", bounding_polygon)
    #print("input", all_points)
    #inpired from: https://stackoverflow.com/questions/36063533/clipping-a-voronoi-diagram-python
    #print("boundary", bounding_polygon)
    #print("to extract", ids_to_extract)
    #print(bounding_geometry)
    #build the voronoi diagram
    #   dt = scipy.spatial.Delaunay(np.array(all_points))
    vor = Voronoi(np.array(all_points), qhull_options ='Qbb Qx')
    centerpoint = np.array(all_points).mean(axis=0)
    #print("Voronoi with x regions and y points", len(vor.regions), len(vor.points))
    #print("input",vor.points)
    #print("vertices", vor.vertices)
    #show what we do
    # from scipy.spatial import voronoi_plot_2d
    # import matplotlib.pyplot as plt
    # fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='black', line_width=3, line_alpha=0.6, point_size=3)
    # plt.show()
    region_outputs = [] #storage for output polygons which will be merged in the last step
    out_stats = []
    #prepare the bounding box geometry
    # border_lines = []
    # for border_pt_idx in range(0, len(bounding_geometry) - 1):
    #     border_lines.append(geometry.LineString([bounding_geometry[border_pt_idx], bounding_geometry[border_pt_idx + 1]]))
    #print("bounding", bounding_polygon)
    #for each of the regions...
    for ids_to_extract in region_dict:
        deviation_dist_list = []
        polygons_list = []  # storage for output polygons which will be merged in the last step of each loop
        reg_stats = []
        #now extract the voronoi cells
        for id in ids_to_extract:
            deviation_dist_list.append(distance_list[id])
            #print("cellpoint", geometry.Point(all_points[id]))
            #print("id", id)
            cell_vertices = vor.regions[vor.point_region[id]] # get the region of the point of which the index is in the input
            #print("cell vertices", cell_vertices)
            polygon_coords = [] #empty list to insert the coordinates of the cell at a later stage
            #border_cell = 0
            added_vertices = [] #store the vertices of the cell for searching ridges in case of border cell
            border_line_intersections = [] #for checking whether the two intersections are on same segment, in case of border cell (if not, need to add corner)
            for v_idx in range(0,len(cell_vertices)):
                if cell_vertices[v_idx] != -1:
                    polygon_coords.append(vor.vertices[cell_vertices[v_idx]])
                    added_vertices.append(v_idx)
                    #print("working on", geometry.Point(vor.vertices[cell_vertices[v_idx]]))
                else:
                    #border_cell = 1
                    if v_idx == 0: #let's just store the previous and next list indicator so we can add far points while keeping the ccw order!
                        prev_v_idx = len(cell_vertices)-1
                        next_v_idx = v_idx+1
                    elif v_idx == len(cell_vertices)-1:
                        prev_v_idx = v_idx-1
                        next_v_idx = 0
                    else:
                        prev_v_idx = v_idx-1
                        next_v_idx = v_idx+1
                    #print(cell_vertices[prev_v_idx], cell_vertices[next_v_idx])
                    pts_to_add = [] #to store the added points and line for retrospective self-intersection check!
                    adding_lines = []
                    if prev_v_idx == next_v_idx:
                        prev_next_vertex = [cell_vertices[prev_v_idx]]
                    else:
                        prev_next_vertex = [cell_vertices[prev_v_idx], cell_vertices[next_v_idx]]
                    #print("prev next", prev_next_vertex)
                    for neighb_vertex in prev_next_vertex: #vertex: using once the previous and once the next vertex
                        for ridge_idx in range(0,len(vor.ridge_vertices)):#to be in ccw order! #scan all the ridge list containing their vertices
                            if vor.ridge_vertices[ridge_idx][0]==-1 and neighb_vertex==vor.ridge_vertices[ridge_idx][1] and id in vor.ridge_points[ridge_idx]: #if one ot them contains the vertex and goes to infinity
                                #print(id)
                                #print("ridges", vor.ridge_vertices[ridge_idx])
                                #print("ridges", vor.ridge_points[ridge_idx])
                                # print("ridges", id, vor.ridge_points[ridge_idx],vor.ridge_vertices[ridge_idx])
                                curr_ridge_pts = vor.ridge_points[ridge_idx] #let's get the two input point indexes that form the line to infinity
                                neighbour_pt = all_points[curr_ridge_pts[1]]
                                #print("neighbours", geometry.Point(neighbour_pt))
                                middle_pt = vor.points[[curr_ridge_pts[0],curr_ridge_pts[1]]].mean(axis=0) #let's calculate a point equidistant to those (order doesn't matter)
                                #print(vor.points[[curr_ridge_pts[0],curr_ridge_pts[1]]], middle_pt)
                                connect_vertex = vor.vertices[vor.ridge_vertices[ridge_idx][1]] #let's use the vertice to which the endless line is connected as a start
                                vector = vor.points[curr_ridge_pts[1]] - vor.points[curr_ridge_pts[0]] #let's calculate a vector perpendicular to the line connecting two input points
                                vector /= np.linalg.norm(vector) #give it a lenght of one
                                vector = np.array([-vector[1], vector[0]])
                                vector = np.sign(np.dot(middle_pt-centerpoint, vector))*vector
                                far_point = connect_vertex + vector*100
                                pts_to_add.append(far_point)
                                adding_lines.append(geometry.LineString([connect_vertex,far_point]))
                                #print("connect", geometry.Point(connect_vertex))
                                #print("far point", far_point)
                                continue
                                #print(geometry.Point(middle_pt))
                                antivector = vector*-1
                                # print("direction", direction)
                                far_point_1 = connect_vertex + vector
                                far_point_2 = connect_vertex + antivector
                                line_1 = geometry.LineString([connect_vertex, far_point_1])
                                line_2 = geometry.LineString([connect_vertex, far_point_2])
                                multiplier = 10
                                while line_1.intersects(bounding_polygon)==False and line_2.intersects(bounding_polygon)==False:
                                    #print("inside while loop")
                                    far_point_1 = far_point_1 + vector * multiplier
                                    far_point_2 = far_point_2 + antivector * multiplier
                                    line_1 = geometry.LineString([connect_vertex, far_point_1])
                                    line_2 = geometry.LineString([connect_vertex, far_point_2])
                                    multiplier = multiplier*10
                                    #print(geometry.Point(connect_vertex), geometry.Point(far_point_1), geometry.Point(far_point_2), line_1, line_2)
                                if line_1.intersects(bounding_polygon)==False:
                                    if bounding_polygon.contains(geometry.Point(connect_vertex))==False:
                                        pts_to_add.append(far_point_1)
                                        adding_lines.append(line_1)
                                        #print("exception caught 1", line_1)
                                    else:
                                        pts_to_add.append(far_point_2)
                                        adding_lines.append(line_2)
                                        #print("2 per def")
                                    continue
                                if line_2.intersects(bounding_polygon)==False:
                                    if bounding_polygon.contains(geometry.Point(connect_vertex))==False:
                                        pts_to_add.append(far_point_2)
                                        adding_lines.append(line_2)
                                        #print("exception caught 2")
                                    else:
                                        pts_to_add.append(far_point_1)
                                        adding_lines.append(line_1)
                                        #print("1 per def")
                                    continue
                                intersect_1 = line_1.intersection(bounding_polygon)
                                #print(intersect_1)
                                dist_1 = (intersect_1.coords[1][0]-connect_vertex[0])**2+(intersect_1.coords[1][1]-connect_vertex[1])**2
                                #print("dist", dist_1, intersect_1)
                                intersect_2 = line_2.intersection(bounding_polygon)
                                #print(intersect_2)
                                dist_2 = (intersect_2.coords[1][0]-connect_vertex[0])**2+(intersect_2.coords[1][1]-connect_vertex[1])**2
                                #print("dist", dist_2, intersect_2)
                                if dist_1>dist_2:
                                    pts_to_add.append(far_point_2)
                                    adding_lines.append(line_2)
                                    #print("2", line_2, line_1, intersect_1, intersect_2)
                                else:
                                    pts_to_add.append(far_point_1)
                                    adding_lines.append(line_1)
                                    #print("1", line_1, line_2, intersect_1, intersect_2, bounding_polygon)
                                #print(far_point)
                    #print("adding lines", adding_lines[0], adding_lines[1])
                    if len(adding_lines)>1 and len(cell_vertices)>2:
                        if adding_lines[0].intersects(adding_lines[1]) is True:
                            print(pts_to_add, adding_lines)
                            new_intersect = adding_lines[0].intersection(adding_lines[1])
                            print("new intersect", new_intersect.coords[0])
                            pts_to_add = [np.array([new_intersect.coords[0][0], new_intersect.coords[0][1]])]
                    #print(pts_to_add)
                    for pt_to_add in pts_to_add:
                        polygon_coords.append(pt_to_add)
                        #print(geometry.Point(pt_to_add))
            #if border_cell==1:
            #print(geometry.Polygon(polygon_coords))
            polygons_list.append(geometry.Polygon(polygon_coords).intersection(bounding_polygon))
            #print(geometry.Polygon(polygon_coords).intersection(bounding_polygon))
            #distance statistics
        mean_distance = np.mean(deviation_dist_list)
        ninety_percentile = np.percentile(deviation_dist_list, 90)
        if ninety_percentile<0:
            ninety_percentile = np.percentile(deviation_dist_list, 10)
        #outstats.append((mean_distance, ninety_percentile))
        out_stats.append((mean_distance, ninety_percentile))
        region_outputs.append(cascaded_union(polygons_list))
    #output = cascaded_union(region_outputs)
    #print("output", output)
    return region_outputs, out_stats
            # if border_cell == 1: # if the cell is a border cell, we need to intersect it with the bounding box which is an input
            #     for ridge_idx in range(0,len(vor.ridge_vertices)): #first, let's find the ridge points of the edges going to infinity
            #         #print("checking", vor.ridge_vertices, vor.ridge_points)
            #         if vor.ridge_vertices[ridge_idx][0]==-1 and id in vor.ridge_points[ridge_idx]:
            #             #print("enter the mess", ridge_idx, vor.ridge_points[ridge_idx])
            #             #print("ridge_vertices", vor.ridge_vertices[ridge_idx])
            #             curr_ridge_pts = vor.ridge_points[ridge_idx]
            #             middle_pt = ((all_points[curr_ridge_pts[1]][0] + all_points[curr_ridge_pts[0]][0]) / 2,
            #                          (all_points[curr_ridge_pts[1]][1] + all_points[curr_ridge_pts[0]][1]) / 2)
            #             equi_pt_1 = (middle_pt[0] - 1000 * (all_points[curr_ridge_pts[1]][1] - all_points[curr_ridge_pts[0]][1]),
            #                          middle_pt[1] + 1000 * (all_points[curr_ridge_pts[1]][0] - all_points[curr_ridge_pts[0]][0]))
            #             equi_pt_2 = (middle_pt[0] + 1000 * (all_points[curr_ridge_pts[1]][1] - all_points[curr_ridge_pts[0]][1]),
            #                          middle_pt[1] - 1000 * (all_points[curr_ridge_pts[1]][0] - all_points[curr_ridge_pts[0]][0]))
            #             # print("points", equi_pt_2, equi_pt_1)
            #             line_voronoi_infinity = geometry.LineString([equi_pt_1, equi_pt_2])
            #             min_dist = np.inf
            #             for border_line_idx in range(0, len(border_lines)):
            #                 intersect = border_lines[border_line_idx].intersection(line_voronoi_infinity)
            #                 #print("intersection", intersect)
            #                 if str(intersect)[0:5]=='POINT':
            #                     sq_dist = (middle_pt[0] - intersect.x) ** 2 + (middle_pt[1] - intersect.y) ** 2
            #                     if sq_dist < min_dist:
            #                         min_dist = sq_dist
            #                         closer_pt = (intersect.x, intersect.y)
            #                         idx_closest_border_line = border_line_idx
            #                         closest_ridge_pts = curr_ridge_pts
            #             #print(idx_closest_border_line)
            #             print("POINT ("+ str(closer_pt[0]) + ' ' + str(closer_pt[1]) + '),')
            #             print(all_points[closest_ridge_pts[1]], all_points[closest_ridge_pts[0]])
            #             polygon_coords.append(closer_pt)
            #             #print("closer pt", closer_pt)
            #             border_line_intersections.append(idx_closest_border_line)
            #     #print(border_line_intersections)
            #     if border_line_intersections[0]!=border_line_intersections[1]:
            #         corner_pts = bounding_geometry[min(border_line_intersections)+1:max(border_line_intersections)+1]
            #         for corner_pt in corner_pts:
            #             polygon_coords.append(corner_pt)
            #print(geometry.MultiPoint(polygon_coords).convex_hull)
                        #print("cell vertices", cell_vertices)
                        # print(geometry.MultiPoint(polygon_coords))
                        # print(geometry.MultiPoint(polygon_coords).convex_hull.intersection(bounding_polygon))
            #polygons_list.append(geometry.MultiPoint(polygon_coords).convex_hull.intersection(bounding_polygon)) #indented
            #print(geometry.Polygon(polygon_coords))
        #print("before intersect", cascaded_union(polygons_list))
        #print(cascaded_union(polygons_list).intersection (bounding_polygon))
        #polygons_list = cascaded_union(polygons_list)#.intersection(bounding_polygon)
        #print(ids_to_extract)
        #print("polygons_list", polygons_list)
    #     region_outputs.append(polygons_list)
    # output = cascaded_union(region_outputs)
    # #output = cascaded_union(region_output)
    # #print(min([item[1] for item in all_points]))
    # #print(min([item[1] for item in vor.vertices]))
    # # try:
    # #     print(max([item[1] for item in output.exterior.coords]))
    # # except:
    # #     print(max([item[1] for item in output.geoms.exterior.coords]))
    # #print("before intersect", output)
    # #print(output)
    # #output = output.intersection(bounding_polygon)
    # print("output", output)
    # return output

# def voronoi_deprecated(all_points, ids_to_extract, bounding_geometry):
    # from shapely.ops import cascaded_union
    # vor = Voronoi(np.array(all_points))
    # #print(vor)
    # # from scipy.spatial import voronoi_plot_2d
    # # import matplotlib.pyplot as plt
    # # fig = voronoi_plot_2d(vor)
    # # plt.show()
    # # print(vor.regions)
    # # print(vor.point_region)
    # # print(vor.vertices)
    # #print(vor.ridge_points)
    # #print(vor.point_region == ids_to_extract)
    # polygons_list = []
    # for id in ids_to_extract:
    #     print("id",id)
    #     print("region_idx", vor.point_region[id])
    #     cell_vertices = vor.regions[vor.point_region[id]]
    #     print("region", cell_vertices)
    #     print("vertices", vor.ridge_vertices)
    #     print("ridges", vor.ridge_points)
    #     polygon_coords = []
    #     pts_to_check = []
    #     for idx in cell_vertices:
    #         print("idx", idx)
    #         if idx==-1:
    #             edges_to_intersect = []
    #             for pt_idx in range(0,len(vor.ridge_points)):
    #                 ridge_v = vor.ridge_vertices[pt_idx]
    #                 ridge_p = vor.ridge_points[idx]
    #                 print(ridge_v, ridge_p)
    #                 #print("pt_idx", pt_idx)
    #                 #print("point index", pt_idx[1])
    #                 if (ridge_p[0]==id or ridge_p[1]==id) and (ridge_[0]==-1 or ridge_v[1]==-1):
    #                     print('enter the mess')
    #                     #print(pt_idx)
    #                     #print("diff", all_points[pt_idx[1]][0]-all_points[pt_idx[0]][0], all_points[pt_idx[1]][1]-all_points[pt_idx[0]][1] )
    #                     middle_pt = ((all_points[pt_idx[1]][0]+all_points[pt_idx[0]][0])/2,(all_points[pt_idx[1]][1]+all_points[pt_idx[0]][1])/2)
    #                     #print("middle_pt", middle_pt)
    #                     equi_pt_1 = (middle_pt[0] -100 * (all_points[pt_idx[1]][1]-all_points[pt_idx[0]][1]), middle_pt[1] +100 * (all_points[pt_idx[1]][0]-all_points[pt_idx[0]][0]))
    #                     equi_pt_2 = (middle_pt[0]+100*(all_points[pt_idx[1]][1]-all_points[pt_idx[0]][1]),middle_pt[1]-100*(all_points[pt_idx[1]][0]-all_points[pt_idx[0]][0]))
    #                     #print("points", equi_pt_2, equi_pt_1)
    #                     line_voronoi_infinity = geometry.LineString([equi_pt_1,equi_pt_2])
    #                     #print(line_voronoi_infinity)
    #                     border_line = geometry.LineString(bounding_geometry)
    #                     #print(border_line)
    #                     intersect = line_voronoi_infinity.intersection(border_line)
    #                     # for border_pt_id in range(0, len(bounding_geometry)-1):
    #                     #     print("border", ([bounding_geometry[border_pt_id],bounding_geometry[border_pt_id+1]]))
    #                     #     border_line = geometry.LineString([bounding_geometry[border_pt_id],bounding_geometry[border_pt_id+1]])
    #                     #     print("intersect", line_voronoi_infinity.intersection(border_line))
    #                     intersect_pts = list(intersect.geoms)
    #                     if len(intersect_pts)>0:
    #                         min_dist = np.inf
    #                         for el in intersect_pts:
    #                             sq_dist = (middle_pt[0]-el.x)**2+(middle_pt[1]-el.y)**2
    #                             if sq_dist<min_dist:
    #                                 min_dist = sq_dist
    #                                 closer_pt = (el.x, el.y)
    #                         #print("pt_to_add", closer_pt)
    #                         polygon_coords.append(closer_pt)
    #                         pts_to_check.append(closer_pt)
    #         else:
    #             #print("vertices", vor.vertices[idx])
    #             polygon_coords.append(vor.vertices[idx])
    #     if len(pts_to_check)==2: #not sur the value can be 1, a voronoi cell on the side should have at least two
    #         intersec_idx_list = [-1, -1]
    #         for pt_check_idx in range(0,len(pts_to_check)):
    #             for border_line_idx in range(0, len(bounding_geometry)-1):
    #                 #print("border", ([bounding_geometry[border_pt_id],bounding_geometry[border_pt_id+1]]))
    #                 border_line = geometry.LineString([bounding_geometry[border_line_idx],bounding_geometry[border_line_idx+1]])
    #                 border_pt_intersec = border_line.intersection(geometry.Point(pts_to_check[pt_check_idx]))
    #                 if str(border_pt_intersec)[0:5]=='POINT':
    #                     intersec_idx_list[pt_check_idx]=border_line_idx
    #         if intersec_idx_list[0]!=intersec_idx_list[1]:
    #             polygon_coords = polygon_coords + bounding_geometry[min(intersec_idx_list)+1:max(intersec_idx_list)+1]
    #             print("indexes", min(intersec_idx_list), max(intersec_idx_list))
    #             print("additional corner", bounding_geometry[min(intersec_idx_list):max(intersec_idx_list)])
    #             print(polygon_coords)
    #     polygons_list.append(geometry.Polygon(polygon_coords))
    # output = cascaded_union(polygons_list)
    # print(output)
    # return output

if __name__ == '__main__':
    coord_list = [[0.86170135,0.58121214],[0.70563184,0.1260651 ],[0.33120399,0.96377443],[0.77696171,0.18824909],[0.7824393,0.8696153 ],[0.84167985,0.23755508],[0.34350368,0.24398097],[0.73284776,0.80344627],[0.94705068,0.26463342],[0.16787806,0.23200608]]
    #np.random.rand(10,2)
    to_extract = [{2: 1, 3: 1, 4: 1}]# 40: 1, 2: 1, 20: 1, 25: 1, 13: 1}
    #coord_list = [[0.5,0.1],[0.1,3.1],[2.1,2.1],[0.1,3.3],[1.1,4.9],[1.5,1.1],[1.1,2.0],[-2.9,3.1],[2.1,0.1],[2.1,3.1],[2.1,9.1],[2.1,4.1]]
    #to_extract = {3:1, 6:1}
    #polygon = [(90769.301, 434678.942), (90766.593, 434679.277), (90766.051, 434674.903), (90768.759, 434674.567),(90769.301, 434678.942)]
    polygon = [(0,0),(1,0),(1,1),(0,1),(0,0)]
    #print(coord_list)
    # points = []
    # for pt in coord_list:
    #     points.append(geometry.Point(pt))
    #voronoi_shape(coord_list, to_extract, polygon, geometry.Polygon(polygon))
    #voronoi_shape([[90630.673, 434544.202], [90630.391, 434544.41500000004], [90630.091, 434544.647], [90629.387, 434545.08400000003], [90631.322, 434544.283], [90631.012, 434544.51300000004], [90630.599, 434544.804], [90630.311, 434545.026], [90630.00600000001, 434545.255], [90632.49100000001, 434543.957], [90632.149, 434544.20900000003], [90631.84, 434544.439], [90631.49, 434544.697], [90631.18000000001, 434544.92600000004], [90630.852, 434545.166], [90630.54400000001, 434545.393], [90630.194, 434545.652], [90633.29000000001, 434543.837], [90632.973, 434544.069], [90632.656, 434544.306], [90632.345, 434544.53500000003], [90632.0, 434544.79000000004], [90631.663, 434545.037], [90631.335, 434545.281], [90630.997, 434545.529], [90630.66900000001, 434545.768], [90630.363, 434545.995], [90630.03, 434546.24], [90634.43400000001, 434543.42100000003], [90634.128, 434543.648], [90633.78, 434543.907], [90633.445, 434544.153], [90633.12, 434544.393], [90632.782, 434544.642], [90632.461, 434544.88], [90632.15000000001, 434545.11], [90631.829, 434545.34500000003], [90631.48700000001, 434545.597], [90631.142, 434545.853], [90630.833, 434546.079], [90630.478, 434546.342], [90630.167, 434546.57], [90636.198, 434542.603], [90635.59300000001, 434543.061], [90635.26, 434543.30700000003], [90634.932, 434543.549], [90634.592, 434543.799], [90634.274, 434544.034], [90633.963, 434544.265], [90633.64, 434544.504], [90633.298, 434544.755], [90632.94900000001, 434545.012], [90632.635, 434545.243], [90632.284, 434545.50200000004], [90631.976, 434545.731], [90631.66, 434545.96400000004], [90631.348, 434546.193], [90630.993, 434546.454], [90630.664, 434546.699], [90630.33, 434546.942], [90637.356, 434542.22000000003], [90637.035, 434542.465], [90636.73700000001, 434542.686], [90636.39, 434542.944], [90636.082, 434543.172], [90635.759, 434543.41000000003], [90635.45300000001, 434543.639], [90635.096, 434543.9], [90634.764, 434544.145], [90634.43400000001, 434544.38800000004], [90634.099, 434544.636], [90633.774, 434544.876], [90633.464, 434545.104], [90633.14600000001, 434545.34], [90632.803, 434545.593], [90632.455, 434545.849], [90632.14600000001, 434546.07800000004], [90631.795, 434546.336], [90631.48700000001, 434546.56200000003], [90631.154, 434546.806], [90638.16100000001, 434542.05700000003], [90637.878, 434542.271], [90637.564, 434542.505], [90637.247, 434542.74], [90636.903, 434542.995], [90636.566, 434543.243], [90636.243, 434543.482], [90635.893, 434543.738], [90635.57800000001, 434543.974], [90635.265, 434544.204], [90634.938, 434544.443], [90634.604, 434544.691], [90634.245, 434544.954], [90633.942, 434545.18], [90633.59300000001, 434545.436], [90633.281, 434545.667], [90632.956, 434545.908], [90632.641, 434546.13800000004], [90632.289, 434546.398], [90639.352, 434541.66000000003], [90639.07400000001, 434541.871], [90638.753, 434542.115], [90638.408, 434542.369], [90638.089, 434542.606], [90637.751, 434542.856], [90637.426, 434543.097], [90637.11200000001, 434543.327], [90636.789, 434543.567], [90636.448, 434543.816], [90636.088, 434544.081], [90635.785, 434544.30700000003], [90635.431, 434544.568], [90635.119, 434544.797], [90634.796, 434545.037], [90634.488, 434545.26300000004], [90634.133, 434545.525], [90633.796, 434545.773], [90633.475, 434546.012], [90633.133, 434546.26300000004], [90640.567, 434541.331], [90640.258, 434541.563], [90639.946, 434541.79600000003], [90639.603, 434542.049], [90639.27500000001, 434542.293], [90638.961, 434542.523], [90638.638, 434542.76300000004], [90638.29400000001, 434543.016], [90637.945, 434543.275], [90637.634, 434543.505], [90637.289, 434543.75800000003], [90636.634, 434544.174], [90636.356, 434544.39], [90636.125, 434544.576], [90635.977, 434544.726], [90635.64, 434544.97500000003], [90635.316, 434545.21400000004], [90634.981, 434545.46], [90634.652, 434545.70300000004], [90641.767, 434540.963], [90641.438, 434541.211], [90641.123, 434541.443], [90640.802, 434541.681], [90640.485, 434541.916], [90640.137, 434542.173], [90639.798, 434542.422], [90639.47200000001, 434542.66500000004], [90639.129, 434542.917], [90638.80500000001, 434543.157], [90638.49100000001, 434543.389], [90638.166, 434543.629], [90637.82400000001, 434543.88], [90637.474, 434544.139], [90637.164, 434544.368], [90636.619, 434544.731], [90636.464, 434544.87700000004], [90636.171, 434545.101], [90635.852, 434545.335], [90642.602, 434540.80100000004], [90642.317, 434541.016], [90641.97200000001, 434541.273], [90641.62, 434541.531], [90641.304, 434541.767], [90640.955, 434542.023], [90640.626, 434542.268], [90640.314, 434542.49700000003], [90639.988, 434542.738], [90639.648, 434542.989], [90639.292, 434543.251], [90638.982, 434543.481], [90638.63100000001, 434543.742], [90638.326, 434543.966], [90637.992, 434544.212], [90637.672, 434544.44800000003], [90637.32400000001, 434544.705], [90636.981, 434544.957], [90636.654, 434545.196], [90643.141, 434540.865], [90642.798, 434541.119], [90642.476, 434541.358], [90642.158, 434541.591], [90641.827, 434541.836], [90641.482, 434542.09], [90641.133, 434542.349], [90640.82, 434542.58], [90640.458, 434542.846], [90640.144, 434543.077], [90639.816, 434543.32], [90639.507, 434543.548], [90639.15000000001, 434543.811], [90638.812, 434544.061], [90638.49100000001, 434544.298], [90638.147, 434544.552], [90637.812, 434544.798], [90644.281, 434540.496], [90643.951, 434540.737], [90643.636, 434540.973], [90643.288, 434541.227], [90642.939, 434541.487], [90642.609, 434541.727], [90642.264, 434541.98600000003], [90641.951, 434542.216], [90641.626, 434542.457], [90641.296, 434542.69800000003], [90640.946, 434542.958], [90640.603, 434543.20900000003], [90640.282, 434543.44800000003], [90639.938, 434543.701], [90639.611, 434543.944], [90639.296, 434544.174], [90638.962, 434544.419], [90646.051, 434539.64], [90645.147, 434540.32300000003], [90644.783, 434540.592], [90644.471, 434540.82300000003], [90644.113, 434541.088], [90643.799, 434541.319], [90643.476, 434541.55700000003], [90643.144, 434541.80100000004], [90642.801, 434542.055], [90642.443, 434542.319], [90642.13100000001, 434542.55], [90641.778, 434542.81], [90641.457, 434543.047], [90641.126, 434543.29000000004], [90640.816, 434543.52], [90640.461, 434543.783], [90640.111, 434544.037], [90639.788, 434544.277], [90646.304, 434539.912], [90645.956, 434540.17], [90645.637, 434540.406], [90645.315, 434540.643], [90644.97200000001, 434540.895], [90644.632, 434541.147], [90644.281, 434541.406], [90643.963, 434541.64], [90643.59, 434541.914], [90643.277, 434542.145], [90642.95300000001, 434542.385], [90642.638, 434542.618], [90642.287, 434542.878], [90641.931, 434543.13800000004], [90641.607, 434543.379], [90641.269, 434543.628], [90640.933, 434543.874], [90648.098, 434539.063], [90647.207, 434539.742], [90646.872, 434539.988], [90646.531, 434540.24], [90646.179, 434540.503], [90645.855, 434540.739], [90645.50200000001, 434541.00200000004], [90645.176, 434541.242], [90644.85, 434541.483], [90644.526, 434541.72000000003], [90644.172, 434541.981], [90643.823, 434542.239], [90643.494, 434542.481], [90643.155, 434542.732], [90642.822, 434542.978], [90642.504, 434543.211], [90642.173, 434543.456], [90641.831, 434543.708], [90649.054, 434538.939], [90648.74, 434539.177], [90648.401, 434539.427], [90648.042, 434539.693], [90647.725, 434539.928], [90647.349, 434540.20300000004], [90647.035, 434540.434], [90646.709, 434540.67600000004], [90646.388, 434540.913], [90646.04000000001, 434541.172], [90645.68400000001, 434541.433], [90645.356, 434541.67600000004], [90645.014, 434541.928], [90644.682, 434542.173], [90644.353, 434542.41500000004], [90644.024, 434542.656], [90643.676, 434542.911], [90643.315, 434543.179], [90650.557, 434538.343], [90649.904, 434538.835], [90649.589, 434539.07], [90649.228, 434539.337], [90648.918, 434539.568], [90648.582, 434539.815], [90648.26, 434540.052], [90647.905, 434540.316], [90647.55, 434540.575], [90646.863, 434541.083], [90646.548, 434541.318], [90646.221, 434541.556], [90645.891, 434541.8], [90645.545, 434542.056], [90645.187, 434542.321], [90644.873, 434542.55100000004], [90644.505, 434542.821], [90651.096, 434538.392], [90650.772, 434538.634], [90650.465, 434538.86], [90650.126, 434539.11], [90649.774, 434539.371], [90649.416, 434539.637], [90649.088, 434539.87700000004], [90648.725, 434540.143], [90648.409, 434540.379], [90648.048, 434540.641], [90647.751, 434540.863], [90647.401, 434541.123], [90647.05500000001, 434541.378], [90646.721, 434541.625], [90646.379, 434541.878], [90646.054, 434542.118], [90645.723, 434542.36100000003], [90645.387, 434542.609], [90652.55900000001, 434537.742], [90651.641, 434538.43700000003], [90651.279, 434538.705], [90650.967, 434538.936], [90650.605, 434539.20300000004], [90650.27500000001, 434539.445], [90649.939, 434539.694], [90649.624, 434539.927], [90649.264, 434540.194], [90648.901, 434540.457], [90648.582, 434540.697], [90648.23300000001, 434540.954], [90647.901, 434541.199], [90647.577, 434541.43700000003], [90647.246, 434541.683], [90646.903, 434541.936], [90646.546, 434542.2], [90646.224, 434542.435], [90652.841, 434538.045], [90652.489, 434538.306], [90652.162, 434538.548], [90651.825, 434538.798], [90651.507, 434539.031], [90651.141, 434539.30100000004], [90650.795, 434539.55700000003], [90650.451, 434539.81], [90650.1, 434540.068], [90649.776, 434540.309], [90649.44900000001, 434540.55], [90649.119, 434540.795], [90648.777, 434541.047], [90648.42, 434541.31200000003], [90648.09300000001, 434541.552], [90647.72200000001, 434541.824], [90647.41, 434542.055], [90654.655, 434537.161], [90654.332, 434537.409], [90654.03, 434537.635], [90653.694, 434537.884], [90653.38, 434538.117], [90653.018, 434538.385], [90652.662, 434538.646], [90652.321, 434538.897], [90651.974, 434539.154], [90651.635, 434539.403], [90651.301, 434539.648], [90650.971, 434539.89400000003], [90650.625, 434540.147], [90650.272, 434540.41000000003], [90649.952, 434540.647], [90649.582, 434540.92], [90649.264, 434541.154], [90648.93400000001, 434541.397], [90648.599, 434541.643], [90655.538, 434536.95900000003], [90655.238, 434537.18700000003], [90654.904, 434537.436], [90654.556, 434537.697], [90654.217, 434537.945], [90653.86200000001, 434538.20900000003], [90653.53600000001, 434538.44800000003], [90653.201, 434538.697], [90652.87700000001, 434538.936], [90652.515, 434539.202], [90652.152, 434539.467], [90651.832, 434539.707], [90651.474, 434539.97000000003], [90651.141, 434540.216], [90650.818, 434540.452], [90650.478, 434540.705], [90650.132, 434540.961], [90649.773, 434541.227], [90649.446, 434541.466], [90656.154, 434537.023], [90655.772, 434537.304], [90655.463, 434537.534], [90655.133, 434537.78], [90654.797, 434538.025], [90654.436, 434538.293], [90654.092, 434538.548], [90653.749, 434538.799], [90653.40000000001, 434539.059], [90653.061, 434539.308], [90652.74100000001, 434539.543], [90652.405, 434539.793], [90652.051, 434540.054], [90651.692, 434540.317], [90651.369, 434540.556], [90651.01, 434540.822], [90650.672, 434541.068], [90657.981, 434536.213], [90657.671, 434536.452], [90657.365, 434536.68200000003], [90657.041, 434536.92600000004], [90656.707, 434537.169], [90656.344, 434537.43700000003], [90655.997, 434537.695], [90655.659, 434537.946], [90655.3, 434538.20900000003], [90654.975, 434538.45], [90654.647, 434538.691], [90654.3, 434538.947], [90653.946, 434539.206], [90653.591, 434539.471], [90653.269, 434539.707], [90652.902, 434539.979], [90652.581, 434540.217], [90652.245, 434540.465], [90651.913, 434540.705], [90658.638, 434536.232], [90658.278, 434536.5], [90657.928, 434536.75800000003], [90657.58, 434537.014], [90657.236, 434537.27], [90656.897, 434537.521], [90656.57, 434537.762], [90656.228, 434538.01300000004], [90655.867, 434538.277], [90655.508, 434538.544], [90655.186, 434538.781], [90654.823, 434539.05100000004], [90654.504, 434539.286], [90654.16500000001, 434539.537], [90653.84, 434539.777], [90653.477, 434540.044], [90653.13, 434540.3], [90659.48300000001, 434536.042], [90659.11200000001, 434536.316], [90658.78600000001, 434536.55700000003], [90658.462, 434536.797], [90658.114, 434537.054], [90657.749, 434537.321], [90657.412, 434537.57300000003], [90657.075, 434537.822], [90656.706, 434538.09], [90656.367, 434538.343], [90656.043, 434538.581], [90655.714, 434538.825], [90655.36600000001, 434539.08400000003], [90654.992, 434539.359], [90654.672, 434539.594], [90654.31, 434539.86100000003], [90653.98300000001, 434540.103], [90660.378, 434535.838], [90660.056, 434536.077], [90659.696, 434536.346], [90659.336, 434536.61], [90658.994, 434536.862], [90658.638, 434537.124], [90658.296, 434537.37700000004], [90657.964, 434537.621], [90657.626, 434537.87], [90657.27500000001, 434538.13], [90656.898, 434538.407], [90656.57800000001, 434538.645], [90656.208, 434538.917], [90655.893, 434539.15], [90655.552, 434539.403], [90655.223, 434539.64400000003], [90654.857, 434539.913], [90661.582, 434535.435], [90661.238, 434535.692], [90660.893, 434535.945], [90660.542, 434536.206], [90660.194, 434536.463], [90659.872, 434536.7], [90659.528, 434536.95300000004], [90659.181, 434537.211], [90658.798, 434537.492], [90658.471, 434537.732], [90658.103, 434538.004], [90657.777, 434538.245], [90657.435, 434538.498], [90657.11, 434538.737], [90656.74100000001, 434539.00800000003], [90656.388, 434539.268], [90656.052, 434539.518], [90662.488, 434535.24], [90662.142, 434535.495], [90661.80900000001, 434535.74], [90661.47, 434535.993], [90661.114, 434536.254], [90660.74, 434536.529], [90660.413, 434536.772], [90660.038, 434537.047], [90659.70300000001, 434537.293], [90659.37, 434537.54000000004], [90659.037, 434537.78500000003], [90658.671, 434538.055], [90658.315, 434538.316], [90657.976, 434538.568], [90657.62, 434538.83], [90657.284, 434539.079], [90656.94, 434539.33], [90663.735, 434534.76900000003], [90663.399, 434535.01900000003], [90663.046, 434535.28], [90662.663, 434535.56200000003], [90662.346, 434535.797], [90661.974, 434536.071], [90661.639, 434536.319], [90661.295, 434536.572], [90660.966, 434536.814], [90660.598, 434537.085], [90660.234, 434537.354], [90659.891, 434537.607], [90659.547, 434537.862], [90659.197, 434538.119], [90658.86600000001, 434538.363], [90658.526, 434538.613], [90658.185, 434538.867], [90665.302, 434534.108], [90664.303, 434534.86], [90663.947, 434535.124], [90663.598, 434535.381], [90663.261, 434535.62700000004], [90662.923, 434535.879], [90662.572, 434536.139], [90662.19, 434536.42], [90661.869, 434536.658], [90661.493, 434536.934], [90661.166, 434537.17600000004], [90660.815, 434537.435], [90660.486, 434537.675], [90660.118, 434537.94800000003], [90659.772, 434538.204], [90659.425, 434538.461], [90666.541, 434533.753], [90665.894, 434534.242], [90665.546, 434534.499], [90665.225, 434534.737], [90664.875, 434534.994], [90664.52100000001, 434535.25800000003], [90664.148, 434535.534], [90663.814, 434535.779], [90663.44900000001, 434536.048], [90662.695, 434536.53], [90663.12700000001, 434536.289], [90662.383, 434536.766], [90662.423, 434536.803], [90662.07, 434537.067], [90661.7, 434537.337], [90661.357, 434537.591], [90661.008, 434537.85000000003], [90660.663, 434538.103], [90667.766, 434533.335], [90667.452, 434533.57300000003], [90666.813, 434534.05100000004], [90666.452, 434534.317], [90666.09300000001, 434534.586], [90665.758, 434534.832], [90665.376, 434535.115], [90665.053, 434535.354], [90664.707, 434535.608], [90663.993, 434536.134], [90663.634, 434536.39900000003], [90662.971, 434536.832], [90662.658, 434537.069], [90662.541, 434537.196], [90662.26, 434537.413], [90661.91, 434537.67100000003], [90661.558, 434537.93], [90668.685, 434533.074], [90668.35, 434533.327], [90668.01, 434533.585], [90667.689, 434533.824], [90667.322, 434534.09500000003], [90666.981, 434534.347], [90666.629, 434534.606], [90666.30900000001, 434534.844], [90665.936, 434535.12], [90665.204, 434535.656], [90664.861, 434535.912], [90664.526, 434536.16000000003], [90664.188, 434536.408], [90663.84, 434536.666], [90663.484, 434536.928], [90663.099, 434537.21], [90662.785, 434537.444], [90662.405, 434537.724], [90669.361, 434532.927], [90669.34300000001, 434532.99700000003], [90669.204, 434533.13800000004], [90668.932, 434533.351], [90668.61200000001, 434533.59], [90668.27100000001, 434533.844], [90667.904, 434534.114], [90667.535, 434534.387], [90667.195, 434534.63800000004], [90666.819, 434534.914], [90666.482, 434535.163], [90666.132, 434535.422], [90665.801, 434535.667], [90665.43000000001, 434535.94], [90665.079, 434536.199], [90664.729, 434536.45900000003], [90664.367, 434536.72500000003], [90664.023, 434536.979], [90663.686, 434537.226], [90669.488, 434533.31], [90669.406, 434533.485], [90669.142, 434533.69], [90668.765, 434533.969], [90668.43400000001, 434534.21400000004], [90668.09, 434534.47000000003], [90667.762, 434534.712], [90667.38100000001, 434534.992], [90667.015, 434535.26], [90666.666, 434535.52], [90666.306, 434535.78500000003], [90665.957, 434536.043], [90665.635, 434536.28], [90665.285, 434536.539], [90664.92, 434536.808], [90664.54000000001, 434537.086], [90669.601, 434533.742], [90669.623, 434533.787], [90669.325, 434534.02], [90668.964, 434534.287], [90668.601, 434534.553], [90668.247, 434534.816], [90667.903, 434535.07], [90667.564, 434535.32], [90667.22200000001, 434535.57300000003], [90666.86600000001, 434535.837], [90666.49100000001, 434536.115], [90666.153, 434536.362], [90665.768, 434536.645], [90669.497, 434534.369], [90669.162, 434534.619], [90668.807, 434534.88300000003], [90668.44, 434535.155], [90668.107, 434535.401], [90667.721, 434535.686], [90667.391, 434535.931], [90667.041, 434536.18700000003], [90666.7, 434536.438], [90669.705, 434534.74], [90669.363, 434534.994], [90669.031, 434535.24], [90668.69, 434535.493], [90668.33, 434535.759], [90667.93400000001, 434536.049], [90669.867, 434535.17], [90669.55500000001, 434535.401], [90669.16500000001, 434535.688], [90630.09, 434547.009], [90630.436, 434546.908], [90630.804, 434546.80700000003], [90631.166, 434546.706], [90631.526, 434546.605], [90631.874, 434546.511], [90632.23, 434546.411], [90632.607, 434546.309], [90632.984, 434546.20300000004], [90633.323, 434546.108], [90633.71800000001, 434546.001], [90634.057, 434545.906], [90634.422, 434545.806], [90634.763, 434545.713], [90635.152, 434545.604], [90635.52100000001, 434545.503], [90635.87, 434545.406], [90636.244, 434545.303], [90636.603, 434545.20300000004], [90636.936, 434545.11100000003], [90637.29400000001, 434545.012], [90637.666, 434544.91000000003], [90638.051, 434544.803], [90638.391, 434544.70900000003], [90638.77100000001, 434544.604], [90639.113, 434544.51], [90639.467, 434544.412], [90639.808, 434544.318], [90629.914, 434546.632], [90630.211, 434546.53500000003], [90630.591, 434546.429], [90630.932, 434546.335], [90631.321, 434546.227], [90631.66500000001, 434546.131], [90632.028, 434546.032], [90632.36200000001, 434545.939], [90632.754, 434545.829], [90633.126, 434545.727], [90633.477, 434545.629], [90633.848, 434545.527], [90634.209, 434545.429], [90634.554, 434545.332], [90634.904, 434545.234], [90635.27100000001, 434545.134], [90635.662, 434545.027], [90636.195, 434544.917], [90636.516, 434544.815], [90636.745, 434544.731], [90637.088, 434544.63300000003], [90637.431, 434544.538], [90637.813, 434544.43200000003], [90638.176, 434544.332], [90638.532, 434544.233], [90638.914, 434544.13], [90639.262, 434544.032], [90639.605, 434543.938], [90639.957, 434543.841], [90640.331, 434543.738], [90640.713, 434543.632], [90641.05500000001, 434543.54000000004], [90641.43000000001, 434543.435], [90641.768, 434543.34], [90642.13100000001, 434543.241], [90642.46, 434543.15], [90642.836, 434543.045], [90643.207, 434542.944], [90643.56, 434542.846], [90643.923, 434542.744], [90644.272, 434542.648], [90644.605, 434542.556], [90644.964, 434542.457], [90645.327, 434542.356], [90645.7, 434542.253], [90646.034, 434542.161], [90646.421, 434542.053], [90646.75600000001, 434541.961], [90647.107, 434541.865], [90647.445, 434541.772], [90647.818, 434541.667], [90648.186, 434541.567], [90648.533, 434541.47000000003], [90648.894, 434541.371], [90649.248, 434541.27400000003], [90649.586, 434541.18], [90649.921, 434541.085], [90650.287, 434540.98600000003], [90650.66500000001, 434540.88], [90651.005, 434540.788], [90651.378, 434540.685], [90651.698, 434540.594], [90652.053, 434540.496], [90652.39, 434540.403], [90652.761, 434540.30100000004], [90653.12700000001, 434540.2], [90653.466, 434540.104], [90653.836, 434540.005], [90654.175, 434539.907], [90654.51, 434539.817], [90654.86200000001, 434539.72000000003], [90655.22, 434539.62200000003], [90655.596, 434539.517], [90655.923, 434539.427], [90656.29400000001, 434539.325], [90656.636, 434539.231], [90656.98, 434539.135], [90657.306, 434539.04600000003], [90657.681, 434538.941], [90658.035, 434538.842], [90658.391, 434538.746], [90658.735, 434538.648], [90659.082, 434538.554], [90659.414, 434538.462], [90659.75600000001, 434538.367], [90660.113, 434538.268], [90660.48700000001, 434538.167], [90660.814, 434538.075], [90661.18400000001, 434537.971], [90661.503, 434537.88300000003], [90661.859, 434537.788], [90662.186, 434537.695], [90662.56, 434537.592], [90662.918, 434537.495], [90663.253, 434537.401], [90663.605, 434537.304], [90663.951, 434537.207], [90664.289, 434537.116], [90664.623, 434537.023], [90664.981, 434536.925], [90665.348, 434536.82300000003], [90665.675, 434536.732], [90666.043, 434536.631], [90666.36200000001, 434536.543], [90666.716, 434536.446], [90667.045, 434536.354], [90667.407, 434536.254], [90667.76, 434536.158], [90668.105, 434536.06200000003], [90668.461, 434535.96400000004], [90668.799, 434535.87200000003], [90669.118, 434535.782], [90629.788, 434546.297], [90629.956, 434546.217], [90630.338, 434546.109], [90630.711, 434546.006], [90631.06, 434545.908], [90631.435, 434545.804], [90631.795, 434545.706], [90632.138, 434545.61100000003], [90632.493, 434545.511], [90632.863, 434545.41000000003], [90633.245, 434545.304], [90633.601, 434545.206], [90633.978, 434545.101], [90634.322, 434545.006], [90634.68000000001, 434544.907], [90635.024, 434544.813], [90635.412, 434544.706], [90636.081, 434544.58], [90636.467, 434544.48], [90636.65000000001, 434544.393], [90636.859, 434544.306], [90637.20300000001, 434544.211], [90637.557, 434544.113], [90637.926, 434544.012], [90638.307, 434543.905], [90638.654, 434543.811], [90639.039, 434543.705], [90639.376, 434543.61100000003], [90639.73700000001, 434543.51300000004], [90640.072, 434543.419], [90640.45, 434543.314], [90640.825, 434543.213], [90641.177, 434543.114], [90641.548, 434543.012], [90641.899, 434542.91500000004], [90642.23, 434542.824], [90642.59300000001, 434542.724], [90642.957, 434542.62200000003], [90643.345, 434542.517], [90643.676, 434542.42600000004], [90644.048, 434542.32], [90644.39, 434542.227], [90644.743, 434542.13], [90645.08, 434542.036], [90645.456, 434541.93200000003], [90645.823, 434541.831], [90646.175, 434541.734], [90646.537, 434541.634], [90646.884, 434541.537], [90647.221, 434541.444], [90647.575, 434541.346], [90647.932, 434541.24700000003], [90648.311, 434541.143], [90648.649, 434541.049], [90649.027, 434540.945], [90649.361, 434540.854], [90649.71, 434540.756], [90650.042, 434540.664], [90650.428, 434540.559], [90650.78, 434540.461], [90651.132, 434540.364], [90651.49, 434540.264], [90651.82800000001, 434540.17], [90652.178, 434540.075], [90652.511, 434539.983], [90652.893, 434539.878], [90653.24100000001, 434539.778], [90653.596, 434539.681], [90653.95300000001, 434539.58400000003], [90654.291, 434539.489], [90654.638, 434539.395], [90654.979, 434539.299], [90655.338, 434539.2], [90655.714, 434539.097], [90656.043, 434539.005], [90656.417, 434538.901], [90656.75, 434538.811], [90657.096, 434538.716], [90657.426, 434538.624], [90657.8, 434538.52], [90658.154, 434538.422], [90658.499, 434538.327], [90658.86, 434538.227], [90659.204, 434538.132], [90659.526, 434538.043], [90659.871, 434537.947], [90660.23300000001, 434537.847], [90660.606, 434537.745], [90660.932, 434537.654], [90661.298, 434537.552], [90661.628, 434537.461], [90661.97, 434537.366], [90662.485, 434537.266], [90662.871, 434537.164], [90663.066, 434537.074], [90663.365, 434536.98], [90663.723, 434536.881], [90664.066, 434536.788], [90664.394, 434536.696], [90664.73300000001, 434536.603], [90665.09, 434536.505], [90665.454, 434536.403], [90665.784, 434536.313], [90666.154, 434536.211], [90666.48, 434536.121], [90666.821, 434536.026], [90667.152, 434535.936], [90667.515, 434535.835], [90667.869, 434535.73600000003], [90668.206, 434535.643], [90668.55900000001, 434535.54600000003], [90668.907, 434535.45], [90669.23700000001, 434535.36100000003], [90669.573, 434535.267], [90669.935, 434535.17], [90629.675, 434545.95900000003], [90629.788, 434545.882], [90630.09700000001, 434545.789], [90630.479, 434545.683], [90630.851, 434545.579], [90631.217, 434545.479], [90631.584, 434545.376], [90631.942, 434545.279], [90632.29000000001, 434545.18200000003], [90632.636, 434545.08400000003], [90633.018, 434544.981], [90633.388, 434544.87700000004], [90633.73, 434544.781], [90634.111, 434544.675], [90634.469, 434544.58], [90634.823, 434544.481], [90635.166, 434544.387], [90635.551, 434544.279], [90635.929, 434544.17600000004], [90636.27500000001, 434544.08], [90636.648, 434543.977], [90637.009, 434543.87700000004], [90637.351, 434543.783], [90637.704, 434543.686], [90638.067, 434543.585], [90638.459, 434543.477], [90638.801, 434543.384], [90639.187, 434543.277], [90639.515, 434543.185], [90639.882, 434543.085], [90640.221, 434542.992], [90640.61, 434542.886], [90640.969, 434542.786], [90641.339, 434542.686], [90641.696, 434542.585], [90642.043, 434542.488], [90642.383, 434542.396], [90642.734, 434542.29600000003], [90643.107, 434542.196], [90643.48700000001, 434542.091], [90643.815, 434541.998], [90644.195, 434541.893], [90644.534, 434541.799], [90644.889, 434541.701], [90645.216, 434541.61100000003], [90645.601, 434541.505], [90645.964, 434541.403], [90646.317, 434541.30700000003], [90646.681, 434541.208], [90647.027, 434541.11], [90647.371, 434541.016], [90647.712, 434540.92], [90648.084, 434540.821], [90648.455, 434540.716], [90648.799, 434540.621], [90649.173, 434540.518], [90649.50200000001, 434540.427], [90649.852, 434540.33], [90650.196, 434540.237], [90650.566, 434540.132], [90650.928, 434540.033], [90651.269, 434539.938], [90651.635, 434539.836], [90651.984, 434539.74], [90652.319, 434539.64900000003], [90652.666, 434539.553], [90653.02100000001, 434539.45300000004], [90653.406, 434539.347], [90653.73, 434539.256], [90654.106, 434539.153], [90654.441, 434539.061], [90654.788, 434538.965], [90655.119, 434538.873], [90655.48700000001, 434538.771], [90655.854, 434538.67100000003], [90656.208, 434538.574], [90656.565, 434538.476], [90656.90000000001, 434538.382], [90657.235, 434538.289], [90657.586, 434538.193], [90657.938, 434538.09500000003], [90658.308, 434537.993], [90658.641, 434537.901], [90659.014, 434537.797], [90659.341, 434537.707], [90659.677, 434537.612], [90660.009, 434537.52], [90660.387, 434537.417], [90660.74, 434537.318], [90661.087, 434537.223], [90661.436, 434537.128], [90661.783, 434537.031], [90662.399, 434536.925], [90662.728, 434536.832], [90663.02, 434536.737], [90663.175, 434536.645], [90663.512, 434536.553], [90663.874, 434536.45300000004], [90664.196, 434536.363], [90664.543, 434536.267], [90664.876, 434536.17600000004], [90665.248, 434536.074], [90665.59300000001, 434535.978], [90665.932, 434535.884], [90666.295, 434535.78500000003], [90666.63, 434535.691], [90666.954, 434535.602], [90667.308, 434535.506], [90667.663, 434535.407], [90668.028, 434535.306], [90668.348, 434535.218], [90668.713, 434535.116], [90669.047, 434535.02400000003], [90669.395, 434534.931], [90669.714, 434534.843], [90629.528, 434545.577], [90629.579, 434545.503], [90629.84700000001, 434545.41500000004], [90630.202, 434545.315], [90630.57800000001, 434545.212], [90630.961, 434545.105], [90631.299, 434545.01], [90631.69, 434544.903], [90632.027, 434544.81], [90632.392, 434544.70900000003], [90632.728, 434544.615], [90633.109, 434544.51], [90633.48700000001, 434544.405], [90633.838, 434544.30700000003], [90634.219, 434544.205], [90634.569, 434544.106], [90634.904, 434544.01300000004], [90635.269, 434543.913], [90635.64600000001, 434543.809], [90636.028, 434543.70300000004], [90636.371, 434543.609], [90636.755, 434543.503], [90637.09300000001, 434543.409], [90637.45, 434543.31], [90637.787, 434543.217], [90638.177, 434543.11], [90638.547, 434543.00800000003], [90638.903, 434542.91000000003], [90639.272, 434542.808], [90639.625, 434542.711], [90639.959, 434542.617], [90640.325, 434542.517], [90640.69, 434542.418], [90641.08, 434542.311], [90641.414, 434542.217], [90641.796, 434542.112], [90642.128, 434542.02], [90642.49100000001, 434541.92100000003], [90642.826, 434541.82800000004], [90643.212, 434541.721], [90643.566, 434541.623], [90643.922, 434541.52400000003], [90644.284, 434541.423], [90644.63100000001, 434541.32800000004], [90644.976, 434541.232], [90645.311, 434541.13800000004], [90645.686, 434541.036], [90646.064, 434540.931], [90646.392, 434540.83900000004], [90646.782, 434540.733], [90647.115, 434540.639], [90647.465, 434540.544], [90647.81, 434540.451], [90648.182, 434540.347], [90648.547, 434540.246], [90648.89, 434540.14900000003], [90649.257, 434540.049], [90649.606, 434539.952], [90649.937, 434539.86100000003], [90650.291, 434539.764], [90650.655, 434539.663], [90651.026, 434539.561], [90651.357, 434539.468], [90651.73, 434539.364], [90652.067, 434539.272], [90652.417, 434539.17600000004], [90652.743, 434539.085], [90653.119, 434538.979], [90653.479, 434538.879], [90654.187, 434538.684], [90654.543, 434538.587], [90654.876, 434538.49700000003], [90655.217, 434538.4], [90655.577, 434538.302], [90655.94900000001, 434538.19800000003], [90656.293, 434538.105], [90656.644, 434538.004], [90656.984, 434537.912], [90657.658, 434537.726], [90658.043, 434537.621], [90658.392, 434537.52400000003], [90658.743, 434537.427], [90659.099, 434537.329], [90659.443, 434537.233], [90659.767, 434537.143], [90660.114, 434537.048], [90660.467, 434536.949], [90660.838, 434536.846], [90661.168, 434536.757], [90661.538, 434536.653], [90661.868, 434536.561], [90662.209, 434536.467], [90662.542, 434536.37700000004], [90662.90000000001, 434536.276], [90663.257, 434536.177], [90663.95300000001, 434535.98600000003], [90664.297, 434535.889], [90664.638, 434535.798], [90664.97, 434535.704], [90665.696, 434535.505], [90666.031, 434535.413], [90666.39600000001, 434535.311], [90666.72, 434535.222], [90667.071, 434535.12700000004], [90667.398, 434535.036], [90667.75600000001, 434534.935], [90668.10800000001, 434534.838], [90668.446, 434534.746], [90668.804, 434534.646], [90669.155, 434534.552], [90669.474, 434534.463], [90669.765, 434534.369], [90629.395, 434545.16500000004], [90629.454, 434545.086], [90629.625, 434545.007], [90630.027, 434544.902], [90630.42, 434544.803], [90630.803, 434544.697], [90631.07800000001, 434544.602], [90631.437, 434544.504], [90631.81, 434544.4], [90632.168, 434544.3], [90632.503, 434544.206], [90632.86200000001, 434544.108], [90633.23, 434544.005], [90633.623, 434543.897], [90633.956, 434543.804], [90634.344, 434543.69800000003], [90634.691, 434543.602], [90635.048, 434543.503], [90635.386, 434543.41000000003], [90635.777, 434543.304], [90636.145, 434543.201], [90636.508, 434543.101], [90636.87, 434543.0], [90637.221, 434542.902], [90637.563, 434542.808], [90637.925, 434542.708], [90638.289, 434542.608], [90638.673, 434542.503], [90639.02100000001, 434542.407], [90639.408, 434542.3], [90639.747, 434542.207], [90640.09700000001, 434542.109], [90640.445, 434542.014], [90640.82800000001, 434541.909], [90641.191, 434541.808], [90641.54400000001, 434541.71], [90641.904, 434541.61100000003], [90642.263, 434541.512], [90642.606, 434541.418], [90642.955, 434541.321], [90643.316, 434541.222], [90643.697, 434541.115], [90644.04000000001, 434541.022], [90644.418, 434540.917], [90644.753, 434540.824], [90645.1, 434540.727], [90645.446, 434540.632], [90645.817, 434540.529], [90646.183, 434540.429], [90646.539, 434540.33], [90646.901, 434540.23], [90647.242, 434540.135], [90647.591, 434540.042], [90647.94, 434539.944], [90648.304, 434539.842], [90648.675, 434539.739], [90649.009, 434539.648], [90649.382, 434539.544], [90649.727, 434539.45], [90650.07400000001, 434539.354], [90650.407, 434539.262], [90650.78600000001, 434539.158], [90651.142, 434539.059], [90651.49, 434538.961], [90651.853, 434538.862], [90652.2, 434538.765], [90652.541, 434538.672], [90652.886, 434538.574], [90653.24, 434538.476], [90653.619, 434538.375], [90653.952, 434538.281], [90654.335, 434538.17600000004], [90654.655, 434538.087], [90655.003, 434537.989], [90655.336, 434537.897], [90655.72200000001, 434537.793], [90656.075, 434537.694], [90656.414, 434537.601], [90656.772, 434537.501], [90657.13100000001, 434537.403], [90657.452, 434537.313], [90657.798, 434537.219], [90658.16, 434537.119], [90658.537, 434537.015], [90658.85800000001, 434536.925], [90659.225, 434536.824], [90659.55900000001, 434536.731], [90659.913, 434536.635], [90660.231, 434536.54600000003], [90660.598, 434536.443], [90660.954, 434536.34500000003], [90661.301, 434536.249], [90661.662, 434536.15], [90662.00200000001, 434536.056], [90662.327, 434535.965], [90662.677, 434535.87], [90663.03, 434535.771], [90663.398, 434535.67], [90663.71800000001, 434535.581], [90664.095, 434535.476], [90664.428, 434535.386], [90664.76, 434535.293], [90665.09, 434535.201], [90665.458, 434535.099], [90665.817, 434535.001], [90666.16, 434534.908], [90666.504, 434534.811], [90666.854, 434534.715], [90667.186, 434534.624], [90667.522, 434534.531], [90667.874, 434534.434], [90668.245, 434534.331], [90668.573, 434534.24], [90668.94, 434534.139], [90669.717, 434534.033], [90669.26, 434534.05100000004], [90629.313, 434544.811], [90629.344, 434544.74], [90629.466, 434544.659], [90629.872, 434544.553], [90630.24, 434544.458], [90630.658, 434544.353], [90630.913, 434544.26], [90635.602, 434542.961], [90635.96800000001, 434542.86], [90636.318, 434542.756], [90636.682, 434542.659], [90637.05, 434542.554], [90637.418, 434542.45900000003], [90637.764, 434542.36100000003], [90638.494, 434542.161], [90638.851, 434542.059], [90639.216, 434541.961], [90639.584, 434541.859], [90639.933, 434541.762], [90640.273, 434541.667], [90640.636, 434541.569], [90641.00600000001, 434541.468], [90641.384, 434541.363], [90641.717, 434541.27], [90642.107, 434541.163], [90642.444, 434541.07], [90642.793, 434540.974], [90643.12700000001, 434540.881], [90643.513, 434540.77400000003], [90643.88100000001, 434540.673], [90644.226, 434540.577], [90644.59300000001, 434540.47500000003], [90644.946, 434540.378], [90645.268, 434540.286], [90645.624, 434540.188], [90645.994, 434540.08900000004], [90646.38100000001, 434539.983], [90646.706, 434539.89], [90647.082, 434539.786], [90647.416, 434539.694], [90647.774, 434539.596], [90648.113, 434539.50200000004], [90648.486, 434539.39900000003], [90648.84700000001, 434539.3], [90649.204, 434539.202], [90649.558, 434539.102], [90649.912, 434539.005], [90650.248, 434538.913], [90650.602, 434538.816], [90650.957, 434538.717], [90651.338, 434538.612], [90651.66100000001, 434538.521], [90652.038, 434538.417], [90652.374, 434538.324], [90652.716, 434538.229], [90653.06, 434538.134], [90653.41900000001, 434538.036], [90653.799, 434537.931], [90654.126, 434537.83900000004], [90654.501, 434537.737], [90654.839, 434537.64400000003], [90655.18400000001, 434537.548], [90655.52500000001, 434537.456], [90655.893, 434537.352], [90656.247, 434537.255], [90656.594, 434537.158], [90656.95300000001, 434537.06], [90657.293, 434536.96400000004], [90657.629, 434536.873], [90657.97, 434536.777], [90658.334, 434536.678], [90658.706, 434536.575], [90659.031, 434536.485], [90659.403, 434536.382], [90659.732, 434536.29000000004], [90660.077, 434536.196], [90660.404, 434536.104], [90660.776, 434536.001], [90661.13100000001, 434535.904], [90661.471, 434535.808], [90661.825, 434535.711], [90662.17, 434535.616], [90662.501, 434535.525], [90662.84300000001, 434535.429], [90663.2, 434535.332], [90663.571, 434535.23], [90663.897, 434535.13800000004], [90664.265, 434535.037], [90664.59300000001, 434534.947], [90664.938, 434534.852], [90665.266, 434534.761], [90665.637, 434534.659], [90665.986, 434534.563], [90666.333, 434534.467], [90666.685, 434534.369], [90667.02500000001, 434534.275], [90667.352, 434534.186], [90667.687, 434534.091], [90668.047, 434533.993], [90668.41, 434533.891], [90668.742, 434533.80100000004], [90669.635, 434533.68], [90638.749, 434541.759], [90639.106, 434541.658], [90639.457, 434541.559], [90639.827, 434541.458], [90640.18000000001, 434541.362], [90640.515, 434541.267], [90640.886, 434541.168], [90641.235, 434541.068], [90641.619, 434540.963], [90641.951, 434540.871], [90642.345, 434540.76300000004], [90642.673, 434540.67100000003], [90646.21800000001, 434539.69], [90647.302, 434539.38800000004], [90647.636, 434539.297], [90647.984, 434539.197], [90648.33, 434539.104], [90651.17, 434538.319], [90651.532, 434538.216], [90653.98700000001, 434537.536], [90654.333, 434537.44], [90654.689, 434537.341], [90655.032, 434537.24700000003], [90655.37700000001, 434537.153], [90655.719, 434537.056], [90656.072, 434536.96], [90656.444, 434536.857], [90657.155, 434536.66000000003], [90657.477, 434536.57], [90657.818, 434536.477], [90658.53, 434536.279], [90658.889, 434536.18200000003], [90659.216, 434536.087], [90659.922, 434535.892], [90660.266, 434535.80100000004], [90660.591, 434535.708], [90660.957, 434535.609], [90661.64, 434535.416], [90662.016, 434535.314], [90662.342, 434535.223], [90662.697, 434535.125], [90663.027, 434535.036], [90663.39, 434534.935], [90663.745, 434534.837], [90664.098, 434534.74], [90664.446, 434534.643], [90665.104, 434534.46], [90665.456, 434534.363], [90665.812, 434534.266], [90666.177, 434534.164], [90666.50200000001, 434534.075], [90666.86600000001, 434533.974], [90667.202, 434533.88], [90667.541, 434533.787], [90667.85800000001, 434533.7], [90668.231, 434533.597], [90668.59, 434533.498], [90669.508, 434533.381], [90668.938, 434533.404], [90663.992, 434534.398], [90664.323, 434534.306], [90665.00200000001, 434534.115], [90665.353, 434534.02], [90665.677, 434533.93], [90666.048, 434533.827], [90666.397, 434533.73], [90666.732, 434533.637], [90667.096, 434533.537], [90667.441, 434533.444], [90667.759, 434533.356], [90668.088, 434533.26], [90668.451, 434533.162], [90669.41900000001, 434533.037]] , [{0: 1, 1265: 1, 1266: 1, 1: 1, 5: 1, 1264: 1, 1149: 1, 1150: 1, 6: 1, 1148: 1, 7: 1, 2: 1, 1147: 1, 1036: 1, 14: 1, 1035: 1, 1037: 1, 13: 1, 8: 1, 1034: 1, 16: 1, 919: 1, 918: 1, 803: 1, 27: 1, 691: 1, 41: 1, 802: 1, 917: 1, 690: 1, 662: 1, 59: 1, 663: 1, 1032: 1, 1144: 1, 3: 1, 1033: 1, 1260: 1, 1261: 1, 1145: 1, 1146: 1, 1262: 1, 1263: 1, 4: 1}, {573: 2, 895: 2, 1010: 2, 894: 2, 1009: 2, 572: 2, 554: 2, 574: 2, 1011: 2, 552: 2, 896: 2, 593: 2, 897: 2, 784: 2, 592: 2, 783: 2, 594: 2, 782: 2, 595: 2, 781: 2, 575: 2, 893: 2, 556: 2, 555: 2, 614: 2, 785: 2, 1012: 2, 786: 2, 613: 2, 787: 2, 612: 2, 901: 2, 788: 2, 630: 2, 789: 2, 629: 2, 790: 2, 791: 2, 643: 2, 792: 2, 793: 2, 652: 2, 794: 2, 906: 2, 641: 2, 907: 2, 651: 2, 908: 2, 795: 2, 796: 2, 797: 2, 658: 2, 798: 2, 657: 2, 799: 2, 800: 2, 801: 2, 661: 2, 913: 2}, {821: 3, 709: 3, 822: 3, 131: 3, 130: 3, 708: 3, 153: 3, 820: 3, 710: 3, 154: 3, 155: 3, 680: 3, 175: 3, 707: 3, 679: 3, 132: 3, 819: 3, 133: 3, 678: 3, 156: 3, 677: 3, 676: 3, 136: 3, 675: 3, 137: 3, 674: 3, 681: 3, 174: 3, 152: 3}, {1258: 4, 1361: 4, 631: 4, 1419: 4, 615: 4, 1434: 4, 653: 4, 644: 4, 1031: 4, 1142: 4, 1141: 4, 645: 4, 1259: 4, 633: 4, 632: 4, 617: 4, 616: 4, 598: 4, 1420: 4, 597: 4, 599: 4, 1360: 4, 600: 4, 1418: 4, 579: 4, 1433: 4, 1417: 4, 578: 4, 1432: 4, 561: 4, 580: 4, 1431: 4, 1416: 4, 581: 4, 1415: 4, 562: 4, 1430: 4, 1429: 4, 1414: 4, 1413: 4, 563: 4, 1412: 4, 543: 4, 1428: 4, 1427: 4, 1426: 4, 1411: 4, 1425: 4, 544: 4, 1410: 4, 1409: 4, 527: 4, 1424: 4, 1408: 4, 1423: 4, 1422: 4, 1421: 4, 1407: 4, 1406: 4, 528: 4, 510: 4, 1405: 4, 1404: 4, 1346: 4, 529: 4, 511: 4, 596: 4}] , [(90670.163, 434535.482), (90629.946, 434547.197), (90629.4810273583, 434545.596432696), (90629.175, 434544.543), (90669.3855851322, 434532.826121509), (90669.684, 434533.845), (90670.163, 434535.482)] , POLYGON Z ((90670.163 434535.482 8.26224097482387, 90629.946 434547.197 8.17874807377439, 90629.4810273583 434545.596432696 7.23630419709833, 90629.175 434544.543 6.61602338048655, 90669.3855851322 434532.826121509 6.69873486617221, 90669.68399999999 434533.845 7.29905731004111, 90670.163 434535.482 8.26224097482387)))
    #print(points)
    #concave_hull, edge_points = alpha_shape(points, alpha=1)
    #print(concave_hull)