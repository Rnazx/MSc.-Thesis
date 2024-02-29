# Set the path as an environment variable. Input the path of the main folder of goblin
export MY_PATH=$(pwd)

#enter the galaxy name
# printf "Choose from the following list of galaxies \n"
# printf "1. M31   ||   2. M33   ||   3. M51   ||   4. NGC6946   \n"
# read galname
# printf "You have selected $galname \n" 

#this makes it case insensitive
galaxy=m33
export galaxy_name=${galaxy,,}

