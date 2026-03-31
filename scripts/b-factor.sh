#! /bin/csh -f

#########################################################
## name your peaklists in Sparky format free.list and bound.list
## name your input pdb file structure.pdb
## the chemical shift perturbations calculated from the two peaklists and normalized to 1H ppm will be incorporated as a b-factor in output.pdb
##########################################################
#creates the chemical shift peturbation file CSP.txt

# paste free.list bound.list | awk '{print $1, (sqrt(((($2-$5)*81.046)^2)+((($3-$6)*799.736)^2)))/799.736}' > _c

# awk '{gsub(/[^0-9. ]/,"")}1' _c > CSP.txt
# rm _c

###########################################################
#creates pdb file with CSP as a b-factor
grep "ATOM" structure.pdb | awk '{printf("%4s %6d %-4s %3s %1s %4d    %7.3f %7.3f %7.3f  %4.2f  %7.6f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,0.0)}'  > structure-0-bfactor.pdb

set res=(`cat CSP.txt  | awk '{print $1}'`)
cat  structure-0-bfactor.pdb > _t.pdb

foreach k ($res)
set bf=(`cat CSP.txt | awk -v nr=$k '{if($1==nr) print $2}'`)
awk -v nres=$k -v bf=$bf '{if($1!="ATOM" || ($1=="ATOM" && $6!=nres)) print $0; if($1=="ATOM" && $6==nres) printf("%4s %6d %-4s %3s %1s %4d    %7.3f %7.3f %7.3f  %4.2f  %7.6f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,bf*800)}' _t.pdb > _tnew.pdb
mv _tnew.pdb _t.pdb
end
awk '{printf("%4s %6d %-4s %3s %1s %4d    %7.3f %7.3f %7.3f  %4.2f  %7.6f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)}' _t.pdb > output.pdb

rm _t.pdb 

echo "created output.pdb, CSP.txt and structure-0-bfactor.pdb" 
