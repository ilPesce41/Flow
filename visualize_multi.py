"""
Script for visualizing flow data in NYC for multiple classes
"""
import geopandas
import geoplot
import pandas as pd
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt
import argparse

def load_data_frame(filename):

    df=pd.read_csv(filename)
    cols_to_keep=['pickup_longitude','pickup_latitude','dropoff_longitude','dropoff_latitude','class']
    df = df[df['pickup_longitude']<-10]
    df = df[df['dropoff_longitude']<-10]
    crs={'init':'epsg:4326'}
    df = df[df.columns.intersection(cols_to_keep)]
    geometry = [LineString([Point(x,y),Point(z,w)]) for x,y,z,w in zip( df["pickup_longitude"],df["pickup_latitude"],df["dropoff_longitude"],df["dropoff_latitude"])]
    df['geometry']=geometry
    geo_df=geopandas.GeoDataFrame(df,crs=crs,)
    return geo_df

def load_single(filename):
    df=pd.read_csv(filename)
    cols_to_keep=['sx','sy','dx','dy','class']

    df = df[df['sy']<-10]
    df = df[df['dy']<-10]
    crs={'init':'epsg:4326'}
    df = df[df.columns.intersection(cols_to_keep)]
    geometry = [LineString([Point(x,y),Point(z,w)]) for x,y,z,w in zip( df["sy"],df["sx"],df["dy"],df["dx"])]
    df['geometry']=geometry
    geo_df=geopandas.GeoDataFrame(df,crs=crs,)
    return geo_df

def main(args):

    boroughs_1 = geopandas.read_file(geoplot.datasets.get_path('nyc_boroughs'))
    boroughs_2 = geopandas.read_file(geoplot.datasets.get_path('nyc_boroughs'))
    all = load_data_frame(args.raw_file)
    single = load_single(args.processed_file)

    fig,ax=plt.subplots(1,2,figsize=(15,15))
    boroughs_1.plot(ax=ax[0],alpha=0.4,color="grey")
    all.plot(ax=ax[0],markersize=200, alpha=0.4,color="green")

    boroughs_2.plot(ax=ax[1],alpha=0.4,color="grey")
    colors = ['red','green','blue','magenta','orange','yello']
    for i,class_val in enumerate(single['class'].unique()):
        single[single['class']==class_val].plot(ax=ax[1],markersize=200, alpha=0.4,color=colors[i])
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_file',type=str)
    parser.add_argument('processed_file',type=str)
    main(parser.parse_args())