import multiprocessing as mp
import gml_pc_merger as utils
import pgSQL_3DCityDB_reader as pg3D
import math

if __name__=='__main__':
    #setting for the number of processes one wishes to run
    process_number = 6
    #get the general data to organize multiprocessing
    inFile = utils.data_loader('globalPC/AHN_buildings.las')
    cityGML_data = pg3D.CityDB_connection("testdb", "127.0.0.1")
    bboxs, ids = cityGML_data.get_all_bboxs()
    intervals = list(range(0, len(ids), math.floor(len(ids)/process_number)))
    intervals.append(len(ids))
    #now that we have divided the points, generate a thread for each group
    for interval_idx in range(0,len(intervals)-1):
        q = mp.Queue()
        start_idx = intervals[interval_idx]
        end_idx = intervals[interval_idx+1]
        process = mp.Process(target=utils.workflow, args=(start_idx, end_idx, bboxs, ids, cityGML_data, inFile, q))
        process.start()
        print(q.get())
