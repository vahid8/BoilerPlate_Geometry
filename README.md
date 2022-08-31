# Content
- [Seraching KDTree](#Serching-KDtree)
- [Cluster points distance base](#Cluster-points-distance-base)
- [Write features to shapefiles](#Write-features-to-shapefiles)
- [Read shapefiles](#Read-shapefiles)
- [Write geojson](#Write-geojson)
- [Read dbf](#Read-dbf)

## Searching KDTree
```
from sklearn.neighbors import KDTree
centers = np.array(centers)
centers_tree = KDTree(centers) # create tree of centers
nearest_dist, nearest_idx = centers_tree.query(cam_6degree["pos"][:2].reshape(1,2), k=4) # search for 4 NEAREST CENTERS to the cam pos
```

## Cluster points distance base
cluster points that are colser than a threshold distance together 
``` 
from sklearn.cluster import DBSCAN
db = DBSCAN(eps=epsilon, min_samples=min_samples).fit(X)
labels = db.labels_
classified_points = {}
for i in range(len(np.unique(labels))):  
    mask = np.where(labels == i)  
    classified_points[i] = X[mask]  
```

## Write features to shapefiles
``` 
from shapely.geometry import Polygon,Point
import geopandas
all_points_3d_left_edge = [[x1,y1,z1], ... , [xn,yn,zn]]
option1:
left_edge_points = list(map(Point, all_points_3d_left_edge))
gdf = geopandas.GeoDataFrame(geometry=left_edge_points)
    if len(gdf) > 0:
        gdf.to_file(os.path.join(output_folder, "shapefiles", "left_edge_points.shp"))
        # or to Geojson file
        # gdf.to_file(os.path.join(output_folder, "shapefiles", "left_edge_points.geojason"), driver='GeoJSON')
option2:
    with fiona.open(os.path.join(output_folder, "shapefiles", "left_edge_line.shp")
            , 'w', 'ESRI Shapefile', {"geometry": "MultiLineString"}) as c:
        ## If there are multiple geometries, put the "for" loop here
        c.write({
            'geometry': mapping(complete_edges),
        })

```


## Read shapefiles
``` 
1) using geopandas
import geopandas
dataframe = geopandas.read_file(test.shp)

2)using shapefile
import shapefile

def read_shape_file(shp_path):
    # -------------------------------------------------
    # Read the shapefiles
    # -------------------------------------------------
    with shapefile.Reader(shp_path) as shp:
        shp_data = list()
        for i in range(shp.numRecords):
            poly_3d = list()
            lines_x_y = shp.shape(i).points
            lines_z = shp.shape(i).z
            bbox = list(shp.shape(i).bbox)
            for num, item in enumerate(lines_x_y):
                poly_3d.append([item[0], item[1], lines_z[num]])

            #Create points in between before add to final result
            poly_3d = expand_points(poly_3d)
            shp_data.append([bbox, np.array(poly_3d, dtype=float)])

        shp_dataframe = pd.DataFrame(shp_data, columns=["bbox", "points"])

    return shp_dataframe

```

## Write geojson 
``` 
import geojson

aa = geojson.Feature(geometry=geojson.Point((100, 100, 0)),properties={"class": "im here"})
bb = geojson.Feature(geometry=geojson.Point((200, 200, 500)),properties={"class": "there"})
feature_collection = geojson.FeatureCollection([aa,bb])

print(feature_collection)
# write to output file
with open("geodata.geojason", 'w') as f:
    geojson.dump(feature_collection, f)
```


#### Read dbf
``` 
from simpledbf import Dbf5
dbf = Dbf5(file_path)
df = dbf.to_dataframe()
```

#### shapely help
####### offset polygons border to the outside
```
coords = [(0, 0),  (0,1), (1, 1), (1, 0) ]
s = Polygon(coords)
t = Polygon(s.buffer(1.0).exterior)
gdf_points = geopandas.GeoDataFrame([s],
                                    columns=["geometry"])
gdf_points.to_file("tt.shp")
gdf_points = geopandas.GeoDataFrame([t],
                                    columns=["geometry"])
gdf_points.to_file("tt2.shp")
```

    

