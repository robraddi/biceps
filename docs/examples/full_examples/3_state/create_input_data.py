import pandas as pd

# Use to reduce large Dataframe to 1 row
#x = pd.read_pickle("3_state/0.noe")
#d = x.to_dict()
#y = [{key:d[key][0] for key in d.keys()}]
#y = pd.DataFrame(y)


# Use to change values in each state
import pandas as pd
state = 0
x = pd.read_pickle("NOE/%s.noe"%(state))
d = x.to_dict()
y = [{key:d[key][0] for key in d.keys()}]
y = pd.DataFrame(y)
y["atom_index1"] = 1
y["atom_name1"] = "H1"
y["atom_index2"] = 20
y["atom_name2"] = "H20"
y["exp"] = 2.5
#model = [3.219, 2.843, 3.133][state]
model = [3.519, 1.4, 3.133][state]
y["model"] = model
y.to_pickle("NOE/%s.noe"%(state))
#y


