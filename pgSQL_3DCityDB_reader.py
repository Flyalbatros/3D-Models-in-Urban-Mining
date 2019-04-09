import psycopg2 as pg

class CityDB_connection(object):

    def __init__(self, dbname, host, password="1234", username="postgres", port="5432"):
        #creates the connexion to the database
        self.conn = pg.connect("dbname={} user={} password={} host={} port={}".format(dbname, username, password, host, port))
        self.cursor = self.conn.cursor()

    def query_and_return(self,query):
        #sends a query to the database and fetches the result
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def get_all_bboxs(self):
        #gets for all buildings in the database the building id and the bounding box ("envelope") as a 3D WKT
        data = self.query_and_return("select c.gmlid, ST_AsText(c.envelope) from cityobject c where c.id in (select b.id from building b where b.building_parent_id is null) union select c.gmlid, ST_AsText(c.envelope) from cityobject c where c.objectclass_id = 25 and c.id not in (select distinct(building_parent_id) from building where building_parent_id is not null);")
        #self.bboxs = [item[0] for item in self.bboxs]
        self.bboxs = []
        self.ids = []
        for el in data:
            self.ids.append(el[0])
            #print(el[1])
            bbox_coords = el[1][10:].replace("(",'').replace(")",'').split(',')
            #original format is wkt string, we convert it here to min max values
            min_x = float(bbox_coords[0].split(' ')[0])
            min_y = float(bbox_coords[0].split(' ')[1])
            #min_z = float(bbox_coords[0].split(' ')[2])
            max_x = float(bbox_coords[1].split(' ')[0])
            max_y = float(bbox_coords[2].split(' ')[1])
            #max_z = float(bbox_coords[2].split(' ')[2])
            self.bboxs.append([(min_x,min_y),(max_x, max_y)])
        return self.bboxs, self.ids

    def get_bbox(self, building_id):
        #gets for all buildings in the database the building id and the bounding box ("envelope") as a 3D WKT
        el = self.query_and_return("SELECT id, ST_AsText(envelope) FROM cityobject where id={} and (objectclass_id=26 or objectclass_id=25);".format(building_id))
        #self.bboxs = [item[0] for item in self.bboxs]
        #print(el[0][1])
        bbox_coords = el[0][1][10:].replace("(",'').replace(")",'').split(',')
        #original format is wkt string, we convert it here to min max values
        min_x = float(bbox_coords[0].split(' ')[0])
        min_y = float(bbox_coords[0].split(' ')[1])
        #min_z = float(bbox_coords[0].split(' ')[2])
        max_x = float(bbox_coords[1].split(' ')[0])
        max_y = float(bbox_coords[2].split(' ')[1])
        #max_z = float(bbox_coords[2].split(' ')[2])
        self.bbox = ([(min_x,min_y),(max_x, max_y)])
        return self.bbox

    def get_roof_geom(self, building_id):
        #for a given buiding, retrieve all the roof geometries (as a 3D WKT)
        roof_geoms = self.query_and_return("SELECT ST_AsText(geometry) from surface_geometry where parent_id in (select lod2_multi_surface_id from citydb_view.thematic_surface where objectclass_id=33 and building_id in (select id from building where building_root_id in (select id from cityobject where gmlid='{}')));".format(building_id))
        #print(roof_geoms)
        #roof_geoms is a list of tuples containing the WKT strings
        out_coord_list = []
        out_wkts = []
        for geom in roof_geoms:
            coord_string_list = geom[0][10:].replace("(",'').replace(")",'').split(',')
            surf_coord_list = []
            for el in coord_string_list:
                coord_string = el.split(' ')
                coord_float = []
                for coord in coord_string:
                    coord_float.append(float(coord))
                surf_coord_list.append(tuple(coord_float))
            out_coord_list.append(tuple(surf_coord_list))
            out_wkts.append(geom[0])
        #print(len(out_coord_list), len(roof_geoms))
        #now we have a list of tuples!
        return out_coord_list, out_wkts

    def get_all_roof_geom(self, id_list):
        self.all_roof_geom = []
        for id in id_list:
            self.all_roof_geom.append(self.get_roof_geom(id))
            # if self.get_roof_geom(id) == []:
            #     print(id)
        return self.all_roof_geom

def wkt_to_coords(wkt_string):
    #print(wkt_string)
    out_coord_list = []
    coord_string_list = wkt_string[wkt_string.index("("):].replace("(", '').replace(")", '').split(', ')
    for el in coord_string_list:
        coord_string = el.split(' ')
        coord_float = []
        for coord in coord_string:
            #print(coord)
            coord_float.append(float(coord))
        out_coord_list.append(tuple(coord_float))
    return out_coord_list

if __name__ == '__main__':
    test = CityDB_connection("testdb", "127.0.0.1")
    bboxs, ids = test.get_all_bboxs()
    all_geom = test.get_all_roof_geom(ids)
    print(bboxs[50])
    #result =  test.query_and_return("SELECT * FROM cityobject LIMIT 10")
    #result = test.get_roof_geom(4834)
    #print(result)
    #for el in result:
        #print(el[0])
