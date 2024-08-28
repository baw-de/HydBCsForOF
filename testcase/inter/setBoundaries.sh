#!/bin/bash
## set -x

echo
echo "##################################################"
echo " Script for setting basic boundary conditions for"
echo " the BAW blueprint for hydraulic structures models"
echo " 2024-08-20  First version  (C. Thorenz)"
echo "##################################################"
echo

NUMBERBCS=0

if [ -e BCS.txt ]  ; then
  rm BCS.txt
fi

while IFS= read -r LINE; do   
    if [ `echo $LINE| grep "(" | wc -l ` -gt 0 ] ; then
      if [ $NUMBERBCS -eq 0 ] ; then
        NUMBERBCS=$OLDLINE
      fi
    fi

    if [ $NUMBERBCS -gt 0 ] ; then
      if [ `echo $LINE| grep "{" | wc -l ` -gt 0 ] ; then
        echo $OLDLINE >> BCS.txt
      fi
    fi

    ## echo $NUMBERBCS
     
    OLDLINE=$LINE
done < constant/polyMesh/boundary 

cp 0/blueprints/head/U 0/bak
cp 0/blueprints/head/alpha.water 0/bak
cp 0/blueprints/head/k 0/bak
cp 0/blueprints/head/omega 0/bak
cp 0/blueprints/head/p_rgh 0/bak
cp 0/blueprints/head/nut 0/bak

for BC in `cat BCS.txt`; do
    ERROR=1
    
    while [ $ERROR -ne 0 ] ; do

      echo 
      echo "###################################################"
      echo "Please choose boundary condition type for \"$BC\":"
      echo 
      echo "   1    Prescribed flowrate"
      echo "   2    Prescribed water level / Air contact"
      echo "   3    Wall (friction)"
      echo "   4    Wall (slip)"
      echo "   5    QUIT"
      read -p "Enter your choice: " BCTYPE
      echo Boundary $BC  Type $BCTYPE  

      if   [  $((BCTYPE)) -eq 1 ] ; then
        BCNAME=flowrate
        read -p "Flowrate [m^3/s]: " FLOWRATE
        ERROR=0
      elif [  $((BCTYPE)) -eq 2 ] ; then  
        BCNAME=waterlevel
        read -p "Water level [m]: " WATERLEVEL
        ERROR=0
      elif [  $((BCTYPE)) -eq 3 ] ; then  
        BCNAME=wall_friction
        read -p "Equivalent sand roughness [mm]: " SANDROUGHNESSMM
        SANDROUGHNESS=`echo $SANDROUGHNESSMM/1000 | bc -l`
        ERROR=0
      elif [  $((BCTYPE)) -eq 4 ] ; then  
        BCNAME=wall_slip
        ERROR=0
      elif [  $((BCTYPE)) -eq 5 ] ; then  
        exit
      else
         echo "###############################################"
         echo "######### ERROR: Wrong number  ################"
         echo "###############################################"
      fi

      if   [  $((ERROR)) -eq 0 ] ; then
        echo "    $BC" >> 0/bak/U 
        echo "    $BC" >> 0/bak/alpha.water
        echo "    $BC" >> 0/bak/k
        echo "    $BC" >> 0/bak/omega 
        echo "    $BC" >> 0/bak/p_rgh
        echo "    $BC" >> 0/bak/nut

        cat 0/blueprints/$BCNAME/k >> 0/bak/k
        cat 0/blueprints/$BCNAME/omega >> 0/bak/omega
        
        cat 0/blueprints/$BCNAME/U.1 >> 0/bak/U
        cat 0/blueprints/$BCNAME/alpha.water.1 >> 0/bak/alpha.water
        cat 0/blueprints/$BCNAME/p_rgh.1 >> 0/bak/p_rgh
        cat 0/blueprints/$BCNAME/nut.1 >> 0/bak/nut

        if   [  $((BCTYPE)) -eq 1 ] ; then
          echo "        flowRateWater $FLOWRATE;        // Constant rate [m³/s]" >> 0/bak/U
        elif [  $((BCTYPE)) -eq 2 ] ; then  
          echo "        data            $WATERLEVEL;                   // water level [m]. Alternative: Supply a table for functions"  >> 0/bak/alpha.water
          echo "        data            $WATERLEVEL;                   // water level [m]. Alternative: Supply a table for functions"  >> 0/bak/p_rgh
        elif [  $((BCTYPE)) -eq 3 ] ; then  
          echo "        Ks              uniform $SANDROUGHNESS;  // roughness height = sand-grain roughness (0 for smooth walls) [m]" >> 0/bak/nut
        fi   

        cat 0/blueprints/$BCNAME/U.2 >> 0/bak/U
        cat 0/blueprints/$BCNAME/alpha.water.2 >> 0/bak/alpha.water
        cat 0/blueprints/$BCNAME/p_rgh.2 >> 0/bak/p_rgh
        cat 0/blueprints/$BCNAME/nut.2 >> 0/bak/nut

        
      fi
   done
done

if   [  $((ERROR)) -eq 0 ] ; then
  echo "}" >> 0/bak/U 
  echo "}" >> 0/bak/alpha.water
  echo "}" >> 0/bak/k
  echo "}" >> 0/bak/omega 
  echo "}" >> 0/bak/p_rgh
  echo "}" >> 0/bak/nut
  cp 0/bak/U 0
  cp 0/bak/alpha.water 0
  cp 0/bak/k 0
  cp 0/bak/omega 0
  cp 0/bak/p_rgh 0
  cp 0/bak/nut 0
 
  echo
  echo "######################################################"
  echo "All files have been copied to the folders 0 and 0/bak."
  echo "Check the boundary parameters in the files in 0/bak." 
  echo "Edit them if necessary and copy the files then "
  echo "to the 0 folder before running setFields!"
  echo "######################################################"
  echo
else
  echo
  echo "###############################################"
  echo "A problem occured. Please try again!"
  echo "###############################################"
  echo
fi

rm BCS.txt





