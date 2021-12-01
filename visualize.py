import geopandas
import geoplot
import pandas as pd
from shapely.geometry import Point, Polygon, LineString
import matplotlib.pyplot as plt

def load_data_frame(filename):
    store_locDF=pd.read_csv(filename)
    cols_to_keep=['pickup_longitude','pickup_latitude','dropoff_longitude','dropoff_latitude']
    store_locDF = store_locDF[store_locDF['pickup_longitude']<-10]
    store_locDF = store_locDF[store_locDF['dropoff_longitude']<-10]
    crs={'init':'epsg:4326'}
    store_locDF = store_locDF[store_locDF.columns.intersection(cols_to_keep)] #Keep only long and lat
    #Convert long and lat to "points"
    geometry = [LineString([Point(x,y),Point(z,w)]) for x,y,z,w in zip( store_locDF["pickup_longitude"],store_locDF["pickup_latitude"],store_locDF["dropoff_longitude"],store_locDF["dropoff_latitude"])]
    #Store geometry in dataframe
    store_locDF['geometry']=geometry
    # geo_df=geopandas.GeoDataFrame(store_locDF,crs=crs,geometry=geopandas.points_from_xy(store_locDF["pickup_longitude"], store_locDF["pickup_latitude"]))
    geo_df=geopandas.GeoDataFrame(store_locDF,crs=crs,)
    return geo_df

def load_single(filename):
    store_locDF=pd.read_csv(filename)
    cols_to_keep=['sx','sy','dx','dy']

    store_locDF = store_locDF[store_locDF['sy']<-10]
    store_locDF = store_locDF[store_locDF['dy']<-10]
    crs={'init':'epsg:4326'}
    store_locDF = store_locDF[store_locDF.columns.intersection(cols_to_keep)] #Keep only long and lat
    #Convert long and lat to "points"
    geometry = [LineString([Point(x,y),Point(z,w)]) for x,y,z,w in zip( store_locDF["sy"],store_locDF["sx"],store_locDF["dy"],store_locDF["dx"])]
    #Store geometry in dataframe
    store_locDF['geometry']=geometry
    # geo_df=geopandas.GeoDataFrame(store_locDF,crs=crs,geometry=geopandas.points_from_xy(store_locDF["pickup_longitude"], store_locDF["pickup_latitude"]))
    geo_df=geopandas.GeoDataFrame(store_locDF,crs=crs,)
    return geo_df

def main():

    boroughs_1 = geopandas.read_file(geoplot.datasets.get_path('nyc_boroughs'))
    boroughs_2 = geopandas.read_file(geoplot.datasets.get_path('nyc_boroughs'))
    all = load_data_frame('tmp1.csv')
    single = load_single("flow_k_out.csv")

    fig,ax=plt.subplots(1,2,figsize=(15,15))
    boroughs_1.plot(ax=ax[0],alpha=0.4,color="grey")
    all.plot(ax=ax[0],markersize=200, alpha=0.4,color="green")

    boroughs_2.plot(ax=ax[1],alpha=0.4,color="grey")
    single.plot(ax=ax[1],markersize=200, alpha=0.4,color="green")
    # plt.legend()
    # fig.show()
    plt.show()
    #https://stackoverflow.com/questions/66498316/plotting-locations-on-a-nyc-map-using-geopandas-and-geoplot

if __name__ == "__main__":
    main()