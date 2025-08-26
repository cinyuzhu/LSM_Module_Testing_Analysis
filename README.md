# LSM_Module_Testing_Analysis
# Analysis scripts meant for LSM module testing 

# Start by creating a folder named Production_Modules_Data
This will host all of the data for high and low temperature

# ------ For analysis of images at high temperature ----- #
cd <Module_data_set>
mkdir Analysis
mkdir Image_1;mkdir Image_2;mkdir Image_3;mkdir Image_4;

# python scripts, these are used for images 1, 2, and 4 mainly for column defects. 
cd Image_1
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_1_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

cd Image_2
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

cd Image_4
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

# ------ For analysis of images at low temperature ----- #
cd <Module_data_set>
mkdir Analysis
mkdir Image_4;

# python scripts, these are used for image 4 mainly for column defects, CTI, noise and sharpness of tracks

cd Image_4
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>






