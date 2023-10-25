#!/usr/bin/env tcsh

if ($#argv == 0 ) then
  echo "usage:    % mfel infile [ outfile [procs [rows [lastp [lastrows] ] ] ] ] " 
  exit 1
endif

set j = 0

echo `hostname`

if ($#argv >= 3) then
  set procs = $3
else
  set procs = `nproc`
endif

if ($#argv >= 4) then
  set rows = $4
else
  set rows = 20
endif

if ($#argv >= 5) then
  set lastp = $5
else
  set lastp = 20
endif

if ($#argv == 6) then
  set lastrows = $6
else
  set lastrows = 5
endif

if($#argv >= 2) then
  echo "mfel $1 $2 $3 $4 $5 $6" >! $2
endif

cp $1 $1.tmpin

while ( $j >= 0  )
   set test1=`grep "^ *eliminate" $1.tmpin`
   set test2=`grep "^ *project" $1.tmpin`
   if("$test1" == "" && "$test2" == "") break

   if($procs <= `nproc`) then
      mpirun -np $procs --oversubscribe -H `hostname` mplrs $1.tmpin -rows rows -lastp lastp -lastrows lastrows $1.tmpout
   else
      mpirun -np $procs --oversubscribe mplrs $1.tmpin -rows rows -lastp lastp -lastrows lastrows $1.tmpout
   endif 

   @ j++
   if($#argv < 2) then
     cat $1.tmpout
   else
     cat $1.tmpout > $2
   endif
   mv -f $1.tmpout $1.tmpin
end
   rm -f $1.tmpin
