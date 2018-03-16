grep "^rule" ~/Documents/new_strain/rules/*/*.rules  | sed "s/.*rules\///"  | sed "s/:/,/" | sed 's/,rule/","/' | sed 's/^/"/' | sed 's/$/"/' |sed "s/://"

