'''
Created in 2025
@author: Liam Thiels

'''
import pandas as pd
import math
import laspy
from laspy.file import File
import numpy as np


#USER INPUT FILES AND PARAMETERS ------------------------------------------------------------------

#This excel needs to have three columns: column 1 "id" is the id of the cross-section, this should have the id for the start and end coordinate. 
#Column 2 "X" and column 3 "Y" needs to contain the X and Y coordinates of the endpoints. So the first row needs to have number 1 in the first column 
#and then the coordinates of the startpoint in column 2 and 3. The second row also needs to have number 1 in the first column and than the coordinates of 
#the endpoint and so forth for all the cross-sections that need to be generated.

Coordinates_path='C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/cross-sectie uiteinden.xlsx' #.xlsx file

#Ideally, the las file has already been clipped to only the area of importance to speed up the processing time.

las_file_path= 'C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/UAVLiDAR_measurements/20241203_Kasterlee_ditches.las' #.las file

Bufferwidth=0.35 #This can still be changed to obtain smaller or wider buffer zones

#Step 1: Load the endpoints of all cross-sections that need to be extracted------------------------

def excel_to_sections(excel_file):
    Endpoints = []
    current_crosssection_id = None
    start_point = None
    
    # Read the Excel file
    df = pd.read_excel(excel_file)

    for index, row in df.iterrows():
        crosssection_id = int(row['id'])  
        x = float(row['X'])  
        y = float(row['Y'])  
        
        if crosssection_id != current_crosssection_id:
            # Start of a new cross-section
            if start_point is not None:
                # If this is not the first cross-section, add the previous cross-section to the list
                Endpoints.append((start_point, end_point))
            start_point = (x, y)
            current_crosssection_id = crosssection_id
        else:
            # Update the endpoint of the cross-section for the current cross-section
            end_point = (x, y)

    if start_point is not None:
        Endpoints.append((start_point, end_point))

    return Endpoints

Endpoints=excel_to_sections(Coordinates_path)
#print(Endpoints)

# Step 2: Calculate the perpendicular distance of all points nearby to the line connecting the endpoints and filter the points that are further than the specified buffer width

# This function calculates the 2D-perpendicular distance of the points in the point cloud to the cross-section
def perpendicular_distance(x,y,transect_start, transect_end):
    x0=x
    y0=y
    # Unpack transect start and end points
    (x1, y1) = transect_start
    (x2, y2) = transect_end
    
    # Calculate the perpendicular distance
    numerator = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))
    denominator = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    
    if denominator == 0:
        # Avoid division by zero in case the transect is a degenerate point
        return float('inf')
    distance= numerator/denominator
    return distance

# This function extracts all positions of points situated within the buffer zone around a cross-section and additionally extracts some metadata of the points
def filter_points(las_file_path, Endpoints, bufferwidth):
    las_file = laspy.read(las_file_path)
    x_coords = las_file['x']
    y_coords = las_file['y']
    z_coords = las_file['z']
    Blue=las_file['blue']
    Green=las_file['green']
    Red=las_file['red']
    Intensity=las_file['intensity']
    GpsTime=las_file['gps_time']
    Angle=las_file['scan_angle_rank']
    Class=las_file['classification']
    filtered_points = []
    for x, y, z, b, g, r, I, gps, a, C in zip(x_coords, y_coords, z_coords, Blue, Green, Red, Intensity, GpsTime, Angle, Class):
        if C ==2: #Classification = 2 --> Ground
         for transect_start, transect_end in Endpoints:
            transect_x_mid = (transect_start[0]+transect_end[0])/2
            transect_y_mid = (transect_start[1]+transect_end[1])/2
            distance_to_mid = np.sqrt((transect_x_mid-x)**2 + (transect_y_mid-y)**2)
            if distance_to_mid < 3:
                distance = perpendicular_distance(x, y, transect_start, transect_end)
                if distance <= Bufferwidth:
                    filtered_points.append((x, y, z, b, g, r, I, gps, a, C))
                    break
    return filtered_points

filtered_points = filter_points(las_file_path, Endpoints, Bufferwidth)
#print(filtered_points)

np.savetxt('C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/UAVLiDAR_measurements/filteredpointsDecember.txt', filtered_points)

# Step 3: Project the points onto the correct cross-section---------------------------------------------------------------------------------------------------------------

#%% Reading the filteredpoints 
filtered_file_path = 'C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/UAVLiDAR_measurements/filteredpointsDecember.txt'

# Initialize an empty list 
filteredpoints = []

# Open the file and process line by line
with open(filtered_file_path, 'r') as file:
    for line in file:
        # Split the line by spaces and convert to floats
        x, y, z, b, g, r, I, gps, a, C = map(float, line.split())
        filteredpoints.append((x, y, z, b, g, r, I, gps, a, C))


#This function projects all the points perpendicularly onto the cross-sectional plane
def project_points_onto_transects(Endpoints, filtered_points):
    projected_points = []

    for transect_start, transect_end in Endpoints:
        for x, y, z, b, g, r, I, gps, a, C in filtered_points:
            transect_x_mid = (transect_start[0] + transect_end[0]) / 2
            transect_y_mid = (transect_start[1] + transect_end[1]) / 2
            distance_to_mid = np.sqrt((transect_x_mid - x) ** 2 + (transect_y_mid - y) ** 2)
            
            if distance_to_mid < 4:
                vector_to_point = np.array([x, y]) - np.array(transect_start)
                line_vector = np.array(transect_end) - np.array(transect_start)
                projection_length = np.dot(vector_to_point, line_vector) / np.dot(line_vector, line_vector)
                if 0<=projection_length<=1:
                    projection = line_vector * projection_length
                    projected_point = np.array(transect_start) + projection
                    projected_points.append((projected_point[0], projected_point[1], z, b, g, r, I, gps, a, C))
    
    return projected_points

projected_points=project_points_onto_transects(Endpoints,filteredpoints)
#print(projected_points)

#Step 4: Assign the points to the correct cross-section------------------------------------------------------------------------------------------------------------------
crossections = {}
for i in range(len(Endpoints)):
    for x, y, z, b, g, r, I, gps, a, C in projected_points:
        transect_x_mid = (Endpoints[i][0][0]+Endpoints[i][1][0])/2
        transect_y_mid = (Endpoints[i][0][1]+Endpoints[i][1][1])/2
        distance_to_mid = np.sqrt((transect_x_mid-x)**2 + (transect_y_mid-y)**2)
        if distance_to_mid < 3: 
            key = f'Transect {i+1}'  # Create the key using string formatting
            if key not in crossections:
                crossections[key] = []  # Initialize the list if it doesn't exist
            crossections[key].append((x, y, z, b, g, r, I, gps, a, C))
#print(crossections)

#%% Save the datapoints with each cross sections
output_lines = []
for transect, points in crossections.items():
    output_lines.append(f'{transect}:')  # Add the transect name
    for point in points:
        output_lines.append(f'{point[0]}, {point[1]}, {point[2]}, {point[3]}, {point[4]}, {point[5]}, {point[6]}, {point[7]}, {point[8]}, {point[9]}') 

# Write the result to a text file
output_path = 'C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/Airborne data/LidarpointspercrosssectioninDecember.txt'
with open(output_path, 'w') as f:
    for line in output_lines:
        f.write(line + '\n')

# Step 5: Convert the 3D coordinates to 2D ------------------------------------------------------------------------------------------------------------------------------

LiDARperXS = 'C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/Airborne data/LidarpointspercrosssectioninDecember.txt'

def load_data(LiDARperXS):
    transects={}
    current_transect=None
    with open(LiDARperXS, 'r')as file:
        for line in file:
            line=line.strip()
            if not line:
                continue
            if line.startswith("Transect"):
                current_transect=line[:-1]
                transects[current_transect]=[]
            elif current_transect:
                try:
                    coordinates=tuple(map(float, line.split(',')))
                    transects[current_transect].append(coordinates)
                except ValueError as e:
                    print(f"Error parsing line'{line}' in {current_transect}:{e}")
    for transect in transects:
        transects[transect]=np.array(transects[transect])
        return  transects
    
crossections=load_data(LiDARperXS)

#This function converts the 3D coordinates to 2D coordinates
def convert_to_2D(crossections):
    transformed_transects = {}
    for transect, points in crossections.items():
        points = np.array(points)
        deepest_idx = np.argmin(points[:, 2]) # Find the index of the deepest point
        deepest_point = points[deepest_idx]
        
        # Compute distance to deepest point in XY plane
        dx = points[:, 0] - deepest_point[0]
        dy = points[:, 1] - deepest_point[1]
        distances = np.sqrt(dx**2 + dy**2)
        distances *= np.sign(dx)  # Left = negative, right = positive

        # Create adjusted points: distance as X, 0 as Y, keep rest unchanged
        adjusted_points = points.copy()
        adjusted_points[:, 0] = distances  
        adjusted_points[:, 1] = 0          

        transformed_transects[transect] = adjusted_points
    return transformed_transects

transformedLiDAR=convert_to_2D(crossections)
print(transformedLiDAR)

def save_all_transects_to_txt(transformed_transects, filename):
    with open(filename, 'w') as f:
        for transect, points in transformed_transects.items():
            f.write(f"{transect}:\n")
            for row in points:
                row_without_y = np.delete(row, 1)  # remove second column (Y = 0)
                row_str = ', '.join(f"{val:.6f}" for val in row_without_y)
                f.write(row_str + '\n')

filename = 'C:/Users/gebruiker/OneDrive - KU Leuven/school/Thesis/Data and scripts/Airborne data/LidarpointspercrosssectioninDecember2D.txt'
save_all_transects_to_txt(transformedLiDAR,filename)
