#!/bin/bash

positions=(182 186 206 208 222 333 337 339 341 351 354 356 445 447 449 482 484 487 488 494)
sequence_string="VPLASRAACEALKDGNGDMVWPNAATVVEVAAWRDAAPATASAAALPEHCEVSGAIAKRTGIDGYPYEIKFRLRMPAEWNGRFFMEGGSGTNGSLSAATGSIGGGQIASALSRNFATIATDGGHDNAVNDNPDALGTVAFGLDPQARLDMGYNSYDQVTQAGKAAVARFYGRAADKSYFIGCSGGREGMMLSQRFPSHYDGIVAGAPGYQLPKAGISGAWTTQSLAPAAVGLDAQGVPLINKSFSDADLHLLSQAILGTCDALDGLADGIVDNYRACQAAFDPATAANPANGQALQCVGAKTADCLSPVQVTAIKRAMAGPVNSAGTPLYNRWAWDAGMSGLSGTTYNQGWRSWWLGSFNSSANNAQRVSGFSARSWLVDFATPPEPMPMTQVAARMMKFDFDIDPLKIWATSGQFTQSSMDWHGATSTDLAAFRDRGGKMILYHGMSDAAFSALDTADYYERLGAAMPGAAGFARLFLVPGMNHCSGGPGTDRFDMLTPLVAWVERGEAPDQISAWSGTPGYFGVAARTRPLCPYPQIARYKGSGDINTEANFACAAPP"
amino_acids=("A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T"
"V" "W" "Y")

# Convert string to array
declare -a sequence=($(echo $sequence_string | grep -o .))

for pindex in "${!positions[@]}"; do
    echo "mutation" >> mutations_${positions[$pindex]}.csv
    for aindex in "${!amino_acids[@]}"; do
        position=${positions[$pindex]}
        new_pindex=$((position - 1))
        echo "${sequence[$new_pindex]}$position${amino_acids[$aindex]}" >> mutations_$position.csv
    done
done

for pindex in "${!positions[@]}"; do
    python esm_dms.py --model-location esm1v_t33_650M_UR90S_1\
    esm1v_t33_650M_UR90S_2 \
    esm1v_t33_650M_UR90S_3 \
    esm1v_t33_650M_UR90S_4 \
    esm1v_t33_650M_UR90S_5 \
    --sequence $sequence_string \
    --dms-input "mutations_${positions[$pindex]}.csv" \
    --dms-output "out_${positions[$pindex]}.csv" \
    --offset-idx 1 \
    --scoring-strategy "masked-marginals" \
    --mutation-col "mutation"
done
