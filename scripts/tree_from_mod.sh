#!/usr/bin/env bash

mod_fn=$1
tree_doctor="/moto/palab/projects/male_mutation_bias_XA/bin/tree_doctor"

# Grab tree string (newick)
tree_str=$( grep "TREE" $mod_fn | awk '{print $NF}' )

# Grab error string and 
err_fn="${mod_fn%.*}.err_est.txt"

if [ -s $err_fn ] ; then
	err_str=$( paste <( $tree_doctor -d $mod_fn | grep dparent | tail -n+2 | awk '{print $NF}' ) <( head -n-12 $err_fn ) | awk -v ORS="," '{scale=$1/$2; print scale * scale * $3}' )
	err_str=${err_str%,} # Drop last comma
else
	err_str="nan"
fi


printf "${tree_str}\t${err_str}\n"
