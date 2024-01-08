# Created by Justin 
# Purpose is to take the files from the dune field calcs, remove the headers, and then change the file names


for file in SpanwiseWSS_*
do grep -v "#" "$file" > temp && mv temp "$file"
done

