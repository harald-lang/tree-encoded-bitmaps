#!/bin/bash

get_mark() {
    local name=$1
    if [[ $name == *"roaring"* ]]; then
         echo "\roaringmark"
    elif [[ $name == *"tree_mask_lo"* ]]; then
         echo "\teblomark"
    elif [[ $name == "tree_mask_po" ]]; then
         echo "\tebpomark"
    elif [[ $name == *"wah"* ]]; then
         echo "\wahmark"
    elif [[ $name == "bitmap" ]]; then
         echo "\uncompressedmark"
    else
        echo "n/a"
    fi
}

get_color() {
    local name=$1
    if [[ $name == *"roaring"* ]]; then
         echo "roaringcolor"
    elif [[ $name == *"tree_mask_lo"* ]]; then
         echo "teblocolor"
    elif [[ $name == "tree_mask_po" ]]; then
         echo "tebpocolor"
    elif [[ $name == *"wah"* ]]; then
         echo "wahcolor"
    elif [[ $name == "bitmap" ]]; then
         echo "uncompressedcolor"
    else
        echo "n/a"
    fi
}

get_print_name() {
    local name=$1
    if [[ $name == *"roaring"* ]]; then
         echo "Roaring"
    elif [[ $name == *"tree_mask_lo"* ]]; then
         echo "TEB-lo"
    elif [[ $name == "tree_mask_po" ]]; then
         echo "TEB-po"
    elif [[ $name == *"wah"* ]]; then
         echo "WAH"
    elif [[ $name == "bitmap" ]]; then
         echo "uncompressed"
    else
        echo "n/a"
    fi
}


DATADIR="./"

SKYLINE_CSVFILE="$DATADIR/skyline_compression_data.csv"
SKYLINE_LOGFILE="$DATADIR/skyline_compression_data.errout"
CSVFILE="$DATADIR/compression_data.csv"
LOGFILE="$DATADIR/compression_data.errout"

#./ex_compression_skyline > $SKYLINE_CSVFILE 2>$SKYLINE_LOGFILE
#./ex_compression > $CSVFILE 2>$LOGFILE

DBFILE="${DATADIR}/compression_data.sqlite3"

#LOAD_MATH="SELECT load_extension('/home/hl/git/sqlite/libsqlitefunctions');"

#export LD_PRELOAD=/lib/x86_64-linux-gnu/libm.so.6
SQL="sqlite3 $DBFILE"

echo "reading $SKYLINE_CSVFILE..."

# a table for the raw results
$SQL "drop table if exists raw_results;"
$SQL "CREATE TABLE raw_results(
    run_id int,
    n int,
    name varchar(20),
    d real,
    d_actual real,
    f real,
    f_actual real,
    size int
    );"
# Import the data
echo ".mode csv
.import $SKYLINE_CSVFILE raw_results" | $SQL


echo "reading $CSVFILE..."

# a table for the raw results
$SQL "drop table if exists raw_results_comp;"
$SQL "CREATE TABLE raw_results_comp(
    run_id int,
    n int,
    name varchar(20),
    d real,
    d_actual real,
    f real,
    f_actual real,
    size int
    );"
# Import the data
echo ".mode csv
.import $CSVFILE raw_results_comp" | $SQL

#exit 0


echo "populating tables with d and f values..."
# Values for d
$SQL "drop table if exists d_values"
$SQL "create table d_values as
        select distinct d as d from raw_results;"
d_s=`$SQL "select d from d_values order by d;"`

$SQL "drop table if exists f_values"
$SQL "create table f_values as
        select distinct f as f from raw_results;"
f_s=`$SQL "select f from f_values order by f;"`

echo -n "values for d are: "
echo  $d_s
echo -n "values for f are: "
echo  $f_s
f_cnt=`echo "$f_s" | wc -l`
echo "considering $f_cnt different values for f."

# adjust y scale for plot
#yscale=`echo "scale=5; 0.0305/$n_cnt" |bc`
xscale="0.7"
yscale="0.7"

# Bitmap types
bitmaps=`$SQL "select distinct name as name from raw_results order by name;"`
echo "bitmaps: $bitmaps"

#--------------------------------------------
# Aggregate the results of the different runs
#--------------------------------------------
echo "aggregating data over the independent runs..."
$SQL "drop table if exists results;"
$SQL "create table results as
  select n, name,
         d, avg(d_actual) as d_actual,
         f, avg(f_actual) as f_actual,
         min(size) size_min,
         avg(size) size_avg,
         max(size) size_max,
         count(*) data_point_cnt
    from raw_results r
   group by n, name, d, f"

$SQL "drop table if exists results_comp;"
$SQL "create table results_comp as
  select n, name,
         d, avg(d_actual) as d_actual,
         f, avg(f_actual) as f_actual,
         min(size) size_min,
         avg(size) size_avg,
         max(size) size_max,
         count(*) data_point_cnt
    from raw_results_comp r
   group by n, name, d, f"


#echo "indexing results..."
#$SQL "create index n_idx on results (n);"

echo "computing skyline..."
$SQL "drop table if exists skyline;"
$SQL "create table skyline as
select d, f,
    json_extract(info, '$.name') as name,
    json_extract(info, '$.size') as size,
    json_extract(info, '$.ratio') as ratio
    from
        (select d_values.d d, f_values.f f, -- cross product
          (select info from
             (select json_object('name', r1.name,
                                 'size', r1.size_avg,
                                 'ratio', (((r1.n  + 7) / 8 + 4) * 1.0) / r1.size_avg
                                 ) as info
                from results r1 where r1.d = d_values.d and r1.f = f_values.f
               order by r1.size_avg limit 1)
          ) as info
         from d_values, f_values
         order by d, f)
;
"

# Create / empty plot file
PLOTFILE="skyline_compression.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0

for f in $bitmaps; do
plotcolor=$(get_color $f)
mark=$(get_mark $f)
tmpfile=`mktemp`
echo ".mode tab
select d, f, 'def' from skyline where name = '$f' order by d, f;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% Bitmap: '$f'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=$mark,xscale=$xscale,yscale=$yscale,$plotcolor}
  }%
]
table[meta=label] {
d f label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

f_name=$(get_print_name $f)
echo "
}; % end of table
\addlegendentry{$f_name}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile

done # for each bitmap type

#exit 0


# Create / empty plot file
PLOTFILE="skyline_compression_ratio.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

f_s=`$SQL "select distinct f from skyline order by f;"`

for f in $f_s; do

echo ".mode tab
select d, f, 1.0/ratio from skyline where name not null and f = $f
union
select d, f, 1.0 from skyline where name is null and f = $f order by f, d;" | $SQL >> $PLOTFILE

echo "" >> $PLOTFILE # empty line
done

#exit 0

# -------------------------------
# Fixed f's
# -------------------------------
f_s='8 64 128'
# Bitmap types
bitmaps=`$SQL "select distinct name as name from raw_results_comp order by name;"`
echo "bitmaps: $bitmaps"
for f in $f_s; do
for b in $bitmaps; do

PLOTFILE="compression_${b}_f${f}.dat"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

echo ".mode tab
select d, size_avg from results_comp where name = '$b' and f = $f order by d;" | $SQL >> $PLOTFILE

done
done




# -------------------------------
# Tree Mask vs WAH
# -------------------------------
echo "computing tree-mask vs WAH skyline..."
$SQL "drop table if exists skyline_treemask_vs_wah;"
$SQL "create table skyline_treemask_vs_wah as
select d, f,
    json_extract(info, '$.name') as name,
    json_extract(info, '$.size') as size,
    json_extract(info, '$.ratio') as ratio
    from
        (select d_values.d d, f_values.f f, -- cross product
          (select info from
             (select json_object('name', r1.name,
                                 'size', r1.size_avg,
                                 'ratio', (((r1.n  + 7) / 8 + 4) * 1.0) / r1.size_avg
                                 ) as info
                from results r1 where r1.d = d_values.d and r1.f = f_values.f
                 and (r1.name like 'tree_mask%' or r1.name like 'wah%' or r1.name like 'bitmap')
               order by r1.size_avg limit 1)
          ) as info
         from d_values, f_values
         order by d, f)
;"

# Create / empty plot file
PLOTFILE="skyline_compression_treemask_vs_wah.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0

for f in $bitmaps; do
plotcolor=$(get_color $f)
mark=$(get_mark $f)
tmpfile=`mktemp`
echo ".mode tab
select d, f, 'def' from skyline_treemask_vs_wah where name = '$f' order by d, f;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% Bitmap: '$f'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=$mark,xscale=$xscale,yscale=$yscale,$plotcolor}
  }%
]
table[meta=label] {
d f label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

f_name=$(get_print_name $f)
echo "
}; % end of table
\addlegendentry{$f_name}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile

done # for each bitmap type

#exit 0


# -------------------------------
# Tree Mask vs Roaring
# -------------------------------
echo "computing tree-mask vs Roaring skyline..."
$SQL "drop table if exists skyline_treemask_vs_roaring;"
$SQL "create table skyline_treemask_vs_roaring as
select d, f,
    json_extract(info, '$.name') as name,
    json_extract(info, '$.size') as size,
    json_extract(info, '$.ratio') as ratio
    from
        (select d_values.d d, f_values.f f, -- cross product
          (select info from
             (select json_object('name', r1.name,
                                 'size', r1.size_avg,
                                 'ratio', (((r1.n  + 7) / 8 + 4) * 1.0) / r1.size_avg
                                 ) as info
                from results r1 where r1.d = d_values.d and r1.f = f_values.f
                 and (r1.name like 'tree_mask%' or r1.name like 'roaring%' or r1.name like 'bitmap')
               order by r1.size_avg limit 1)
          ) as info
         from d_values, f_values
         order by d, f)
;"

# Create / empty plot file
PLOTFILE="skyline_compression_treemask_vs_roaring.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"
plotcolor=$(get_color $f)
mark=$(get_mark $f)
for f in $bitmaps; do
plotcolor=$((plotcolor + 1))

tmpfile=`mktemp`
echo ".mode tab
select d, f, 'def' from skyline_treemask_vs_roaring where name = '$f' order by d, f;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% Bitmap: '$f'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=$mark,xscale=$xscale,yscale=$yscale,$plotcolor}
  }%
]
table[meta=label] {
d f label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

f_name=$(get_print_name $f)
echo "
}; % end of table
\addlegendentry{$f_name}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile

done # for each bitmap type

#exit 0


# -------------------------------
# Roaring vs WAH
# -------------------------------
echo "computing tree-mask vs Roaring skyline..."
$SQL "drop table if exists skyline_roaring_vs_wah;"
$SQL "create table skyline_roaring_vs_wah as
select d, f,
    json_extract(info, '$.name') as name,
    json_extract(info, '$.size') as size,
    json_extract(info, '$.ratio') as ratio
    from
        (select d_values.d d, f_values.f f, -- cross product
          (select info from
             (select json_object('name', r1.name,
                                 'size', r1.size_avg,
                                 'ratio', (((r1.n  + 7) / 8 + 4) * 1.0) / r1.size_avg
                                 ) as info
                from results r1 where r1.d = d_values.d and r1.f = f_values.f
                 and (r1.name like 'wah%' or r1.name like 'roaring%' or r1.name like 'bitmap')
               order by r1.size_avg limit 1)
          ) as info
         from d_values, f_values
         order by d, f)
;"

# Create / empty plot file
PLOTFILE="skyline_compression_roaring_vs_wah.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0

for f in $bitmaps; do
plotcolor=$(get_color $f)
mark=$(get_mark $f)
tmpfile=`mktemp`
echo ".mode tab
select d, f, 'def' from skyline_roaring_vs_wah where name = '$f' order by d, f;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% Bitmap: '$f'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=$mark,xscale=$xscale,yscale=$yscale,$plotcolor}
  }%
]
table[meta=label] {
d f label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

f_name=$(get_print_name $f)
echo "
}; % end of table
\addlegendentry{$f_name}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile

done # for each bitmap type

exit 0



echo "extracting blocked bloom filter results..."
$SQL "drop table if exists bbf;"
$SQL "create table bbf as
select
  json_extract(filter, '$.word_size') * json_extract(filter, '$.w') as block_size,
  json_extract(filter, '$.word_size') as word_size,
  json_extract(filter, '$.w') as word_cnt,
  json_extract(filter, '$.s') as sector_cnt,
  json_extract(filter, '$.z') as zone_cnt,
  (json_extract(filter, '$.word_size') * json_extract(filter, '$.w')) / json_extract(filter, '$.s') as sector_size,
  json_extract(filter, '$.k') as k,
  json_extract(filter, '$.addr') as addr,
  case when json_extract(filter, '$.w') <= json_extract(filter, '$.s')
    then
        case when json_extract(filter, '$.w') = 1
        then 'single'
        else
         case when json_extract(filter, '$.z') > 1
         then 'zoned'
         else 'seq'
         end
        end
    else 'rnd'
  end as access,
  m, b, n, false_positives, fpr, cycles_per_lookup, filter
 from results
where json_extract(filter, '$.name') = 'blocked_bloom_multiword';
"
$SQL "create index n_idx_bbf on bbf (n);"


echo "computing blocked bloom filter skyline..."
$SQL "drop table if exists bbf_skyline;"
$SQL "create table bbf_skyline as
select n_values.n, tw_values.tw,
  (select filter_info from
    (select json_object('block_size',bbf.block_size,
                        'sector_size',bbf.sector_size,
                        'word_size',bbf.word_size,
                        'w', bbf.word_cnt,
                        's', bbf.sector_cnt,
                        'z', bbf.zone_cnt,
                        'k', bbf.k,
                        'access', bbf.access,
                        'addr', bbf.addr,
                        'overhead', (bbf.cycles_per_lookup + bbf.fpr * tw_values.tw),
                        'size', json_extract(bbf.filter, '$.size')
                        ) as filter_info,
            (bbf.cycles_per_lookup + bbf.fpr * tw_values.tw) as overhead
       from bbf
      where bbf.n = n_values.n
      order by overhead limit 1
    )
  ) as filter_info
 from n_values, tw_values;"


echo "determining block sizes..."
bs=`$SQL "select distinct json_extract(filter_info, '$.block_size') as block_size from bbf_skyline order by block_size;"`
echo "block sizes:"
echo $bs

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_blocksize.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0
for b in $bs; do
plotcolor=$((plotcolor + 1))

tmpfile=`mktemp`
echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.block_size') = $b order by tw, n;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% block size: '$b'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,blocksizecolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$b Bytes\,\,}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile
done


# Block Access Patterns

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_access_pattern.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

#a_s=`$SQL "select distinct json_extract(filter_info, '$.access') as info from bbf_skyline order by info;"`
a_s="rnd
seq
zoned
single"

# do not distinguish between sequential and zoned
#a_s=`$SQL "select distinct case when json_extract(filter_info, '$.access') = 'zoned' then 'seq' else json_extract(filter_info, '$.access')  end as info from bbf_skyline order by info;"`
echo "access patterns are: $a_s"

plotcolor=0
for a in $a_s; do
plotcolor=$((plotcolor + 1))
echo "
% access: '$a'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,accesscolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.access') = '$a' order by tw, n;" | $SQL >> $PLOTFILE
#echo ".mode tab
#select tw, n, 'def' from bbf_skyline where case when json_extract(filter_info, '$.access') = 'zoned' then 'seq' else json_extract(filter_info, '$.access') end = '$a' order by tw, n;" | $SQL >> $PLOTFILE

a_name="n/a"
if [[ $a == "seq" ]]; then
#    a_name="Sequential"
    a_name="Sectorized"
elif [[ $a == "rnd" ]]; then
    a_name="Blocked"
#    a_name="Random"
elif [[ $a == "single" ]]; then
    a_name="Register-blocked"
#    a_name="Single word"
elif [[ $a == "zoned" ]]; then
    a_name="Cache-sectorized"
#    a_name="Sequential/Zoned"
fi
echo "
}; % end of table
\addlegendentry{$a_name}
" | sed 's/_/ /g' >> $PLOTFILE

done


# Sector count
# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_sector_cnt.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

sc=`$SQL "select distinct json_extract(filter_info, '$.s') as sector_cnt from bbf_skyline order by sector_cnt;"`

plotcolor=0
for s in $sc; do
plotcolor=$((plotcolor + 1))
echo "
% sector cnt: '$s'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,sectorcntcolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.s') = $s order by tw, n;" | $SQL >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$s\,\,}
" | sed 's/_/ /g' >> $PLOTFILE

done


# Zone count
# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_zone_cnt.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

zc=`$SQL "select distinct json_extract(filter_info, '$.z') as zone_cnt from bbf_skyline order by zone_cnt;"`

plotcolor=0
for z in $zc; do
plotcolor=$((plotcolor + 1))
echo "
% zone cnt: '$z'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,zonecntcolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.z') = $z order by tw, n;" | $SQL >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{\$z=$z\$\,\,}
" | sed 's/_/ /g' >> $PLOTFILE

done


# Word size
# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_word_size.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

ws=`$SQL "select distinct json_extract(filter_info, '$.word_size') as word_size from bbf_skyline order by word_size;"`

plotcolor=0
for w in $ws; do
plotcolor=$((plotcolor + 1))
echo "
% word size: '$w'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,wordcntcolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.word_size') = $w order by tw, n;" | $SQL >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$w Bytes\,\,}
" | sed 's/_/ /g' >> $PLOTFILE

done


# K's

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_k.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

k_s=`$SQL "select distinct json_extract(filter_info, '$.k') as info from bbf_skyline order by info;"`

for k in $k_s; do
echo "
% k: '$k'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,kcolor$k}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.k') = $k order by tw, n;" | $SQL >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$k}
" >> $PLOTFILE

done


# Addressing modes

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_addr.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

#addr_s=`$SQL "select distinct json_extract(filter_info, '$.addr') as info from bbf_skyline order by info;"`

plotcolor=0
#for addr in $addr_s; do
for addr in pow2 magic; do
plotcolor=$((plotcolor + 1))
echo "
% addr: '$addr'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,addrcolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.addr') = '$addr' order by tw, n;" | $SQL >> $PLOTFILE

addr_name="n/a"
if [[ $addr == "pow2" ]]; then
    addr_name="Power of two"
elif [[ $addr == "magic" ]]; then
    addr_name="Magic"
fi
echo "
}; % end of table
\addlegendentry{$addr_name}
" >> $PLOTFILE

done




# -------------------------------
# Cuckoo Skyline
# -------------------------------

echo "extracting cuckoo filter results..."
$SQL "drop table if exists cf;"
$SQL "drop index if exists n_idx_cf;"
$SQL "create table cf as
select
  json_extract(filter, '$.tag_bits') as tag_bits,
  json_extract(filter, '$.associativity') as associativity,
  json_extract(filter, '$.addr') as addr,
  m, b, n, false_positives, fpr,cycles_per_lookup, filter
 from results
where json_extract(filter, '$.name') = 'cuckoo';
"
echo "indexing cuckoo filter results..."
$SQL "create index n_idx_cf on cf (n);"


echo "computing cuckoo filter skyline..."
$SQL "drop table if exists cf_skyline;"
$SQL "create table cf_skyline as
select n_values.n, tw_values.tw,
  (select filter_info from
    (select json_object('tag_bits', cf.tag_bits,
                        'associativity', cf.associativity,
                        'overhead', (cf.cycles_per_lookup + cf.fpr * tw_values.tw),
                        'size', json_extract(filter, '$.size'),
                        'addr', cf.addr
                        ) as filter_info,
            (cf.cycles_per_lookup + cf.fpr * tw_values.tw) as overhead
       from cf
      where cf.n = n_values.n
      order by overhead limit 1
    )
  ) as filter_info
 from n_values, tw_values;"
# -------------------------------


# -------------------------------
# Compute the performance differences (BBF vs CF)
# -------------------------------
echo "computing relative speedups (cuckoo vs. blocked bloom)..."
$SQL "drop table if exists perf_diff;"
$SQL "create table perf_diff as
select bs.n as n, bs.tw as tw,
  case
    when json_extract(bs.filter_info, '$.overhead') < json_extract(cs.filter_info, '$.overhead')
    then 'bbf'
    else 'cf'
  end as filter_class,
  case
    when json_extract(bs.filter_info, '$.overhead') < json_extract(cs.filter_info, '$.overhead')
    then json_extract(cs.filter_info, '$.overhead') / json_extract(bs.filter_info, '$.overhead')
    else json_extract(bs.filter_info, '$.overhead') / json_extract(cs.filter_info, '$.overhead')
  end as speedup
   from bbf_skyline bs, cf_skyline cs
  where bs.n=cs.n and bs.tw=cs.tw
  order by n, tw;
"

PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_bbf_cf_performance_diff.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

inf=99999999
declare -a speedup_steps=(1 1.05 1.1 1.25 1.5 1.75 2 3 4 5 10 $inf)

len=${#speedup_steps[@]}

for ((i=1; i<${len}; i++)); do

begin_idx=$((i-1))
end_idx=i
begin=${speedup_steps[$begin_idx]}
end=${speedup_steps[$end_idx]}

echo "
% s: '[$begin, $end)'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,speedupcolor$i}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from perf_diff where speedup >= $begin and speedup < $end order by tw, n;" | $SQL >> $PLOTFILE

label="[$begin, $end)"
if [[ $end == $inf ]]; then
    label="\$\ge\$ $begin\,x"
elif [[ $end == 1.05 ]]; then
    label="\$<\$ 5\,\%"
fi
echo "
}; % end of table
\addlegendentry{$label}
" >> $PLOTFILE

done
# -------------------------------


# -------------------------------
# Cuckoo tag size
ts=`$SQL "select distinct json_extract(filter_info, '$.tag_bits') as tag_bits from cf_skyline order by json_extract(filter_info, '$.tag_bits');"`

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_cf_tag_size.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0
for t in $ts; do
plotcolor=$((plotcolor + 1))

tmpfile=`mktemp`
echo ".mode tab
select tw, n, 'def' from cf_skyline where json_extract(filter_info, '$.tag_bits') = $t order by tw, n;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% tag size: '$t'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,blocksizecolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$t Bits\,\,}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile
done
# -------------------------------


# -------------------------------
# Cuckoo associativity
as=`$SQL "select distinct json_extract(filter_info, '$.associativity') as associativity from cf_skyline order by associativity;"`

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_cf_associativity.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0
for a in $as; do
plotcolor=$((plotcolor + 1))

tmpfile=`mktemp`
echo ".mode tab
select tw, n, 'def' from cf_skyline where json_extract(filter_info, '$.associativity') = $a order by tw, n;" | $SQL > $tmpfile
outputsize=`cat $tmpfile | wc -l`
if [[ $outputsize != 0 ]]; then
echo "
% associativity: '$a'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,blocksizecolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

cat $tmpfile >> $PLOTFILE

echo "
}; % end of table
\addlegendentry{$a \,\,}
" | sed 's/_/ /g' >> $PLOTFILE
fi
rm $tmpfile
done

# -------------------------------




# -------------------------------
# Memory footprint
# -------------------------------

PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_memory_footprint.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

inf=1073741824
declare -a sizes=(  0       32768     1048576    10485760         67108864         134217728         268435456            $inf)
declare -a labels=("" "$\le$\,L1" "$\le$\,L2" "$\le$\,L3" "$\le$\,64\,MiB" "$\le$\,128\,MiB" "$\le$\,256\,MiB" "$>$\,256\,MiB")

len=${#sizes[@]}

for ((i=1; i<${len}; i++)); do

begin_idx=$((i-1))
end_idx=i
begin=${sizes[$begin_idx]}
end=${sizes[$end_idx]}

echo "
% size: '[$begin, $end)'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,memfootprintcolor$i}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from bbf_skyline where json_extract(filter_info, '$.size') > $begin and json_extract(filter_info, '$.size') <= $end order by tw, n;" | $SQL >> $PLOTFILE

label=${labels[$i]}
echo "
}; % end of table
\addlegendentry{$label}
" >> $PLOTFILE

done
# -------------------------------
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_cf_memory_footprint.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

inf=1073741824
#declare -a sizes=(  0       32768      262144    36700160         67108864         134217728         268435456            $inf)
declare -a sizes=(  0       32768     1048576    10485760         67108864         134217728         268435456            $inf)
declare -a labels=("" "$\le$\,L1" "$\le$\,L2" "$\le$\,L3" "$\le$\,64\,MiB" "$\le$\,128\,MiB" "$\le$\,256\,MiB" "$>$\,256\,MiB")

len=${#sizes[@]}

for ((i=1; i<${len}; i++)); do

begin_idx=$((i-1))
end_idx=i
begin=${sizes[$begin_idx]}
end=${sizes[$end_idx]}

echo "
% size: '[$begin, $end)'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,memfootprintcolor$i}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from cf_skyline where json_extract(filter_info, '$.size') > $begin and json_extract(filter_info, '$.size') <= $end order by tw, n;" | $SQL >> $PLOTFILE

label=${labels[$i]}
echo "
}; % end of table
\addlegendentry{$label}
" >> $PLOTFILE

done
# -------------------------------


# -------------------------------
# Addressing modes
# -------------------------------

# Create / empty plot file
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_cf_addr.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

plotcolor=0
#for addr in $addr_s; do
for addr in pow2 magic; do
plotcolor=$((plotcolor + 1))
echo "
% addr: '$addr'
\addplot[
  scatter,
  only marks,
  point meta=explicit symbolic,
  scatter/classes={
    def={mark=square*,xscale=$xscale,yscale=$yscale,addrcolor$plotcolor}
  }%
]
table[meta=label] {
tw n label
" >> $PLOTFILE

echo ".mode tab
select tw, n, 'def' from cf_skyline where json_extract(filter_info, '$.addr') = '$addr' order by tw, n;" | $SQL >> $PLOTFILE

addr_name="n/a"
if [[ $addr == "pow2" ]]; then
    addr_name="Power of two"
elif [[ $addr == "magic" ]]; then
    addr_name="Magic"
fi
echo "
}; % end of table
\addlegendentry{$addr_name}
" >> $PLOTFILE

done
# -------------------------------



# -------------------------------
# Memory footprint for selected n's
# -------------------------------

PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_memory_footprint_detail.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

#declare -a sizes=(       115097       920781     10417458)
#declare -a labels=("\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$")

declare -a sizes=(       115097       920781     10417458    103496016)
declare -a labels=("\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$" "\$n=10^8\$")

#declare -a sizes=(        12098        96785       920781     12388515)
#declare -a labels=("\$n=10^4\$" "\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$")


#declare -a sizes=(         9741       115097       920781     10417458    103496016)
#declare -a labels=("\$n=10^4\$" "\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$" "\$n=10^8\$")

len=${#sizes[@]}

for ((i=0; i<${len}; i++)); do

n=${sizes[$i]}

color_idx=$((i+1))
echo "
% n: '$n'
\addplot[
  draw=plotcolor$color_idx
]
table  {
tw size
" >> $PLOTFILE

echo ".mode tab
select tw, json_extract(filter_info, '$.size')/1024.0/1024.0 from bbf_skyline where n = $n order by tw;" | $SQL >> $PLOTFILE

label=${labels[$i]}
echo "
}; % end of table
\addlegendentry{$label}
" >> $PLOTFILE

done
# -------------------------------
PLOTFILE="/home/hl/git/papers/bloomfilter/skyline/skyline_cf_memory_footprint_detail.tex"
echo -n "" > $PLOTFILE
echo "writing plot file: `basename $PLOTFILE`"

#declare -a sizes=(       115097       920781     10417458)
#declare -a labels=("\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$")
declare -a sizes=(       115097       920781     10417458    103496016)
declare -a labels=("\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$" "\$n=10^8\$")
#declare -a sizes=(         9741       115097       920781     10417458    103496016)
#declare -a labels=("\$n=10^4\$" "\$n=10^5\$" "\$n=10^6\$" "\$n=10^7\$" "\$n=10^8\$")

len=${#sizes[@]}

for ((i=0; i<${len}; i++)); do

n=${sizes[$i]}

color_idx=$((i+1))
echo "
% n: '$n'
\addplot[
  draw=plotcolor$color_idx
]
table  {
tw size
" >> $PLOTFILE

echo ".mode tab
select tw, json_extract(filter_info, '$.size')/1024.0/1024.0 from cf_skyline where n = $n order by tw;" | $SQL >> $PLOTFILE

label=${labels[$i]}
echo "
}; % end of table
\addlegendentry{$label}
" >> $PLOTFILE

done
# -------------------------------






#rm $DBFILE
echo "Database file: $DBFILE"

# --- Impala vs ours ---
# Our implementation (w=s=k=8) performs not as good as the Impala implementation.
# Reason: The latter does not not use gather instructions.
#
#drop table tcost;
#drop table t;
#
#-- pre-compute the overhead (costs)
#create table tcost as select r.*, tw, (r.cycles_per_lookup + r.fpr * tw_values.tw) as overhead from results r, n_values, tw_values where r.n = n_values.n;
#
#filter                                       filter_name           m           n           b           false_positives   fpr         num_data_points  lookups_per_sec  cycles_per_lookup  tw          overhead
#-------------------------------------------  --------------------  ----------  ----------  ----------  ----------------  ----------  ---------------  ---------------  -----------------  ----------  ----------
#{"name":"blocked_bloom_impala","size":1024}  blocked_bloom_impala  8192        1024        8           2912444.11764706  0.0325428   17               275825000.0      8.69111            8           8.9514524
#{"name":"blocked_bloom_impala","size":1024}  blocked_bloom_impala  8192        1024        8           2912444.11764706  0.0325428   17               275825000.0      8.69111            16          9.2117948
#{"name":"blocked_bloom_impala","size":1024}  blocked_bloom_impala  8192        1024        8           2912444.11764706  0.0325428   17               275825000.0      8.69111            32          9.7324796
#{"name":"blocked_bloom_impala","size":1024}  blocked_bloom_impala  8192        1024        8           2912444.11764706  0.0325428   17               275825000.0      8.69111            64          10.7738492
#{"name":"blocked_bloom_impala","size":1024}  blocked_bloom_impala  8192        1024        8           2912444.11764706  0.0325428   17               275825000.0      8.69111            128         12.8565884
#
#
#-- find out where Impala is faster
#create table t as select distinct n, tw from skyline where json_extract(filter,'$.name') like '%impala';
#
#
#-- find the minimum overhead for each filter and n,tw
#drop table best_per_filter; create table best_per_filter as select n, tw, json_extract(filter,'$.name') as filter_name, min(overhead) as overhead, min(cycles_per_lookup) cycles_per_lookup from tcost group by n, tw, json_extract(filter,'$.name') order by n, tw;
#
#
#-- compute the delta
#select t.n, t.tw, impala.overhead as impala_overhead, bbf.overhead as bbf_overhead, (1.0 - impala.overhead / bbf.overhead) * 100 as percent, bbf.overhead - impala.overhead as cycles from t, best_per_filter impala, best_per_filter bbf where impala.filter_name = 'blocked_bloom_impala' and t.n = impala.n and t.tw = impala.tw and bbf.filter_name = 'blocked_bloom_multiword' and t.n = bbf.n and t.tw = bbf.tw order by cycles desc;
#n           tw          impala_overhead   bbf_overhead      percent           cycles
#----------  ----------  ----------------  ----------------  ----------------  ----------------
#262144      65536       15.0975492047059  17.0213457317647  11.3022586896211  1.92379652705882
#524288      65536       15.3687947482353  17.2870032        11.0962462930805  1.91820845176471
#4194304     65536       14.8437942682353  16.7269927247059  11.2584401001687  1.88319845647059
#1048576     65536       15.0921812988235  16.9481857223529  10.951050772836   1.85600442352941
#2097152     65536       14.96929536       16.78603632       10.8229299959003  1.81674096
#131072      65536       14.31706384       15.93737184       10.1667201861559  1.620308
#4194304     32768       13.7061471341176  14.9480565741176  8.30816657564937  1.24190944
#2097152     32768       13.79904768       14.9108341552941  7.45623258709088  1.11178647529412
#65536       65536       12.5833549929412  13.63893024       7.7394284484501   1.05557524705882
#131072      131072      16.69802768       17.6039269835294  5.1460069357081   0.90589930352941
#262144      131072      17.4931984094118  18.3211933647059  4.51932872936747  0.82799495529412
#524288      131072      17.7541894964706  18.57621712       4.4251615827874   0.82202762352941
#4194304     131072      17.1190885364706  17.9076191152941  4.40332449415388  0.78853057882352
#1048576     131072      17.4284625976471  18.2011238117647  4.2451291585535   0.77266121411764
#1048576     32768       13.9240406494118  14.6588672752941  5.01284725540029  0.73482662588235
#2097152     131072      17.30979072       17.97039504       3.67607010602476  0.66060431999999
#4194304     16384       13.1373235670588  13.3989782870588  1.95279605947802  0.26165472000000
#524288      32768       14.1760973741176  14.3386738823529  1.13383224675599  0.16257650823529
#2097152     16384       13.21392384       13.3066670776471  0.69696819726445  0.09274323764705
#262144      32768       13.8997246023529  13.92920752       0.21166256303329  0.02948291764705
#32768       65536       11.08802856       11.1006129176471  0.11336633157484  0.01258435764706


# --- Sequential vs Random block access pattern ---
#




# ---------
# Early out
# ---------

# !!! requires a special executable !!!


# SIMD calibration disabled to prevent unrolling (which reduces the probability of an early out)
#for s in 0.001 0.01 0.1 1.0; do                                 echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked32 MULTI_WORD_CNT_LO=2 MULTI_WORD_CNT_HI=16 MULTI_SECTOR_CNT_LO=2 MULTI_SECTOR_CNT_HI=16 K_LO=2 K_HI=16 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.0025 0.025 0.25 0.005 0.05 0.5 0.0075 0.075 0.75; do echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked32 MULTI_WORD_CNT_LO=2 MULTI_WORD_CNT_HI=16 MULTI_SECTOR_CNT_LO=2 MULTI_SECTOR_CNT_HI=16 K_LO=2 K_HI=16 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.001 0.01 0.1 1.0; do                                 echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked64 MULTI_WORD_CNT_LO=2 MULTI_WORD_CNT_HI=8  MULTI_SECTOR_CNT_LO=2 MULTI_SECTOR_CNT_HI=8  K_LO=2 K_HI=16 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.0025 0.025 0.25 0.005 0.05 0.5 0.0075 0.075 0.75; do echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked64 MULTI_WORD_CNT_LO=2 MULTI_WORD_CNT_HI=8  MULTI_SECTOR_CNT_LO=2 MULTI_SECTOR_CNT_HI=8  K_LO=2 K_HI=16 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done

## TODO parse raw file and import
#
## a table for the raw results
#$SQL "CREATE TABLE early_out_raw_results(
#    filter json,
#    m int,
#    b int,
#    n int,
#    s real,
#    insert_time_nanos int,
#    false_positives int,
#    fpr real,
#    lookups_per_sec real,
#    cycles_per_lookup real,
#    thread_cnt int,
#    scalar_code int
#    );"
#
#CSVFILE="$DATADIR/early_out.csv"
#
#echo ".mode csv
#.import $CSVFILE early_out_raw_results
#select count(*) from early_out_raw_results;" | $SQL
#
#
#$SQL "create view early_out_s_values as select distinct s from early_out_raw_results;"
#
##--------------------------------------------
## Aggregate the results of the different runs
##--------------------------------------------
## Performance results
## Notes:
##  - scalar_code is deprecated
##  - performance measurements can be identified with n=0
##  - WARNING: remove the unrolling attribute u before group by filter!
#$SQL "drop view early_out_results;
#create view early_out_results as
#  select json_remove(filter,'$.u', '$.e') as filter,
#         json_extract(filter, '$.name') as filter_name,
#         m, n, b, s, thread_cnt,
#         avg(false_positives) as false_positives,
#         avg(fpr) as fpr,
#         max(lookups_per_sec) as lookups_per_sec,
#         min(cycles_per_lookup) as cycles_per_lookup,
#         count(*) as num_data_points
#    from early_out_raw_results r
#   where b<=32 --and fpr>0
#   group by json_remove(filter,'$.u'), m, n, b, s, thread_cnt;"
#
#$SQL "drop view early_out_bbf;
#create view early_out_bbf as
#select
#  json_extract(filter, '$.word_size') * json_extract(filter, '$.w') as block_size,
#  json_extract(filter, '$.word_size') as word_size,
#  json_extract(filter, '$.w') as word_cnt,
#  json_extract(filter, '$.s') as sector_cnt,
#  (json_extract(filter, '$.word_size') * json_extract(filter, '$.w')) / json_extract(filter, '$.s') as sector_size,
#  json_extract(filter, '$.k') as k,
#  json_extract(filter, '$.addr') as addr,
#  case when json_extract(filter, '$.w') <= json_extract(filter, '$.s')
#    then
#        case when json_extract(filter, '$.w') = 1
#            then 'single'
#            else 'seq'
#        end
#    else 'rnd'
#  end as access,
#  m, b, n, s, cycles_per_lookup, filter
# from early_out_results
#where json_extract(filter, '$.name') = 'blocked_bloom_multiword';
#"
#
#$SQL "create view early_out_n_values as select distinct n from early_out_results;"
#
#$SQL "drop view early_out_bbf_skyline;
#create view early_out_bbf_skyline as
#select early_out_s_values.s, early_out_n_values.n, tw_values.tw,
#  (select filter_info from
#    (select json_object('block_size', eobbf.block_size,
#                        'sector_size', eobbf.sector_size,
#                        'word_size', eobbf.word_size,
#                        'w', eobbf.word_cnt,
#                        's', eobbf.sector_cnt,
#                        'k', eobbf.k,
#                        'access', eobbf.access,
#                        'addr', eobbf.addr,
#                        'overhead', (eobbf.cycles_per_lookup + bbf.fpr * tw_values.tw)) as filter_info, -- take the fpr from bbf table to reduce noise
#            (eobbf.cycles_per_lookup + bbf.fpr * tw_values.tw) as overhead
#       from early_out_bbf eobbf, bbf
#      where eobbf.s = early_out_s_values.s
#        and eobbf.filter=bbf.filter
#        and eobbf.n=bbf.n
#        and eobbf.m=bbf.m
#        and bbf.n=early_out_n_values.n
#      order by overhead limit 1
#    )
#  ) as filter_info
# from early_out_n_values, early_out_s_values, tw_values
# where tw_values.tw<=131072;
# select * from early_out_bbf_skyline;
# "
#
## should be empty ( empirical evidence => early out in combination with multiword blocking is slower than branch free)
#$SQL " select json_extract(s.filter_info, '$.overhead') - json_extract(eos.filter_info, '$.overhead') as delta, * from early_out_bbf_skyline eos, bbf_skyline s where eos.n = s.n and eos.tw = s.tw and ( json_extract(s.filter_info, '$.overhead') - json_extract(eos.filter_info, '$.overhead')) > 0 order by delta;"


# Early out in combination with register-blocking
# doesn't make sense, as k <= 4 ???
#for s in 0.001 0.01 0.1 1.0; do                                 echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked32 MULTI_WORD_CNT_LO=1 MULTI_WORD_CNT_HI=1 MULTI_SECTOR_CNT_LO=1 MULTI_SECTOR_CNT_HI=1 K_LO=2 K_HI=8 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.0025 0.025 0.25 0.005 0.05 0.5 0.0075 0.075 0.75; do echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked32 MULTI_WORD_CNT_LO=1 MULTI_WORD_CNT_HI=1 MULTI_SECTOR_CNT_LO=1 MULTI_SECTOR_CNT_HI=1 K_LO=2 K_HI=8 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.001 0.01 0.1 1.0; do                                 echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked64 MULTI_WORD_CNT_LO=1 MULTI_WORD_CNT_HI=1 MULTI_SECTOR_CNT_LO=1 MULTI_SECTOR_CNT_HI=1 K_LO=2 K_HI=8 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done
#for s in 0.0025 0.025 0.25 0.005 0.05 0.5 0.0075 0.075 0.75; do echo "Running with s=$s"; BENCH_PERFORMANCE=1 BENCH_PRECISION=0 RUNS=1 FILTERS=multiregblocked64 MULTI_WORD_CNT_LO=1 MULTI_WORD_CNT_HI=1 MULTI_SECTOR_CNT_LO=1 MULTI_SECTOR_CNT_HI=1 K_LO=2 K_HI=8 THREAD_CNT_LO=1 THREAD_CNT_HI=1 THREAD_STEP_MODE=1 THREAD_STEP=1 BITS_PER_ELEMENT_HI=32 N_LO=1048576 N_HI=1048576 SIMD_CALIBRATION=0 SEL=$s numactl -N 0 ./skyline_earlyout 2>> skyline_earlyout_$s.errout >> skyline_earlyout_$s.out; done

#sqlite> select min(k_per_word), avg(k_per_word), max(k_per_word) from (select (json_extract(filter_info, '$.k') * 1.0) / json_extract(filter_info, '$.w') as k_per_word from bbf_skyline order by k_per_word asc);
#min(k_per_word)  avg(k_per_word)   max(k_per_word)
#---------------  ----------------  ---------------
#0.25             1.58000469924812  6.0