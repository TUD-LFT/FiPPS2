#!/bin/bash

declare -a arr=("input_naca2412" "input_naca0012" "input_naca4415" "input_nlf0115" "input_fx77w121" "input_e423" "input_clarky"  "input_airfoil_optExample" "input_airfoil_optExample_2")


# declare -a arr=("input_airfoil_optExample_2" "input_airfoil_optExample" "input_clarky" "input_e423" "input_fx77w121" "input_naca0012" "input_naca0030" "input_naca2412" "input_naca4210" "input_naca4415" "input_nlf0115" "input_nlf1015")

# declare -a arr=("input_nlf0115")

rm "run_all_ca-cw.dat"
touch "run_all_ca-cw.dat"
for i in "${arr[@]}"; do
    echo $i >> "run_all_ca-cw.dat"
    for j in $(seq -2 10); do
        t=`echo $j | sed -e 's#\.#,#g'`
        printf -v s "%23.16E" $t
        s=`echo $s | sed -e 's#,#\.#g'`
        sed -i "6s# 1.0000000000000000E+99#"$s"#g" ""$i".dat"
        ../PANEL2D $i
        sed -i "6s#"$s"# 1.0000000000000000E+99#g" ""$i".dat"
        alfa=`grep "aoa" $i'_results.dat' | colrm 25`
        ca=`grep "C_L" $i'_results.dat' | colrm 25`
        cw=`grep "C_D" $i'_results.dat' | colrm 25`
        echo $alfa" "$ca" "$cw >> "run_all_ca-cw.dat"
    done
    for j in $(seq 15 15); do
        t=`echo $j | sed -e 's#\.#,#g'`
        printf -v s "%23.16E" $t
        s=`echo $s | sed -e 's#,#\.#g'`
        sed -i "6s# 1.0000000000000000E+99#"$s"#g" ""$i".dat"
        ../PANEL2D $i
        sed -i "6s#"$s"# 1.0000000000000000E+99#g" ""$i".dat"
        alfa=`grep "aoa" $i'_results.dat' | colrm 25`
        ca=`grep "C_L" $i'_results.dat' | colrm 25`
        cw=`grep "C_D" $i'_results.dat' | colrm 25`
        echo $alfa" "$ca" "$cw >> "run_all_ca-cw.dat"
    done
    echo " " >> "run_all_ca-cw.dat"
    rm $i"_*.dat"
done
