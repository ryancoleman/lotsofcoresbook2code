a=41
for i in s*; do
  new=$(printf "s%d" "$a") 
  cp -R -- "$i" "$new"
  let a=a+1
done
