#!/bin/bash

# Prepare file

pip install -r requirements.txt

# Prepare directory
rm -rf ../data/*


# Run file
echo "This script will run the following program to simulate ranked gene lists and RSM scores"
python3 main.py -h


# Ask if the user wants to use default parameters 
read -p "Do you want to use default parameters?(Y/N) " usedefault 
if [ $usedefault == "Y" ]; then
	python3 main.py 
else
	read -p "Please input the options following the command python3 main.py " options
	python3 main.py $options
fi

echo
echo "The program is completed! Thank you"
# python3 main.py
