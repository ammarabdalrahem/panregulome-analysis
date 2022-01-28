#!usr/bin/env bash


#Identif the variables

# add your username and password
username="login=mohamedabdelfadeel1994@gmail.com"
password="password=Rido@1994"
file=$1 #?
#file="Bdistachyon_283_v2.1.gene.gff3.gz"
organism="BrachyPan"

#Login

curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode $username --data-urlencode $password -c cookies > /dev/null


#Retrieve all

curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=BrachyPan' -b cookies > $organism.file_list

#cat $organism.file_list

#Filter names
data=$(grep $file $organism.file_list | cut -d ' ' -f 13 | cut -b 5- | sed "s/amp;//" | cut -c 3- | rev | cut -c2- |rev)
	
echo $data

#Get data


curl  'https://genome.jgi.doe.gov/'$data -b cookies > $file













