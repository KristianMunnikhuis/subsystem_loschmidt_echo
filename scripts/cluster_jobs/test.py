my_list = [[1, 2, 3],[3,4,5]]

# write each item on its own line
with open("list.txt", "w") as f:
    for item in my_list:
        f.write(f"{item}\n")

# or save the entire list as JSON
import json
with open("list.json", "w") as f:
    json.dump(my_list, f)
