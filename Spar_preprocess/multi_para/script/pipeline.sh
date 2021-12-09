
declare -a S_others=("CBS432" "N44" "UFRJ50816" "YPS138" "UWOPS919171")


for S_other in "${S_others[@]}"
do 
	./align_to_Scer.sh $S_other
	./lastz_to_chain_net_axt.sh $S_other
done
