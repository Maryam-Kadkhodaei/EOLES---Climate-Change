"""
30 year test_add_climate

Import modules and libraries needed for the programme
"""
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd


# We set the name of the model here, it will be used in outputs name
model_name = ""

"""INITIALISATION OF THE MODEL"""

model = pyo.ConcreteModel()

# Dual Variable, used to get the marginal value of an equation.
model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

"""INPUTS"""
time_step = 2
num_of_years = 30

# Production profiles of VRE
#load_factor = pd.read_csv("inputs_30y/vre_{}_{}_resampled_2h.csv".format(model_name, scenario), index_col=[0, 1], header=None)
#load_factor = load_factor.squeeze().copy()

# Demand profile in each our in GW
#demand = pd.read_csv("inputs_30y/demand_historical_30y_2h.csv", index_col=0, header=None)
#demand = demand.squeeze().copy()

# Monthly lake inflows in GWh
lake_inflows = pd.read_csv("inputs_30y/lake_inflow_30years.csv", index_col=0, header=None)
lake_inflows = lake_inflows.squeeze().copy()

# Existing capacities of the technologies by December 2017 in GW
capa_ex = pd.read_csv("inputs_30y/existing_capas_elec_new.csv", index_col=0, header=None)
capa_ex = capa_ex.squeeze().copy()

# Maximum capacities of the technologies in GW
capa_max = pd.read_csv("inputs_30y/max_capas_elec_new.csv", index_col=0, header=None)
capa_max = capa_max.squeeze().copy()

# Fixed capacities of the technologies in GW
fix_capa = pd.read_csv("inputs_30y/fix_capas.csv", index_col=0, header=None)
fix_capa = fix_capa.squeeze().copy()

# Annualized power capex cost in M€/GW/year
capex = pd.read_csv("inputs_30y/annuities_elec_new.csv", index_col=0, header=None)
capex = capex.squeeze().copy()

# Annualized energy capex cost of storage technologies in M€/GWh/year
capex_en = pd.read_csv("inputs_30y/str_annuities_elec_new.csv", index_col=0, header=None)
capex_en = capex_en.squeeze().copy()

# Annualized fixed operation and maintenance costs M€/GW/year
fOM = pd.read_csv("inputs_30y/fO&M_elec_new.csv", index_col=0, header=None)
fOM = fOM.squeeze().copy()

# Variable operation and maintenance costs in M€/GWh
vOM = pd.read_csv("inputs_30y/vO&M_elec_new.csv", index_col=0, header=None)
vOM = vOM.squeeze().copy()

# Charging related annuity of storage in M€/GW/year
s_capex = pd.read_csv("inputs_30y/s_capex.csv", index_col=0, header=None)
s_capex = s_capex.squeeze().copy()

# Charging related fOM of storage in M€/GW/year
s_opex = pd.read_csv("inputs_30y/s_opex.csv", index_col=0, header=None)
s_opex = s_opex.squeeze().copy()

# Charging efficiency of storage technologies
eta_in = pd.read_csv("inputs_30y/eta_in.csv", index_col=0, header=None)
eta_in = eta_in.squeeze().copy()

# Discharging efficiency of storage technolgoies
eta_out = pd.read_csv("inputs_30y/eta_out.csv", index_col=0, header=None)
eta_out = eta_out.squeeze().copy()

# Existing storage capacity in GWh
capacity_ex = pd.read_csv("inputs_30y/capacity_ex.csv", index_col=0, header=None)
capacity_ex = capacity_ex.squeeze().copy()
# Parameters in miscellaneous.csv :
# eta_ccgt                  : efficiency of CCGT power plants with CCS
# H2_demand                 : hourly hydrogen demand on top of the storage
# eta_electrolysis          : efficiency of Electrolysis
# phs_discharging_lower     : lower bounds for capa(phs)
# phs_discharging_upper     : upper bounds for capa(phs)
# phs_charging_lower        : lower bounds for s(phs)
# phs_charging_upper        : upper bounds for s(phs)
# phs_energy_lower          : lSower bounds for capacity(phs)
# phs_energy_upper          : upper bounds for capacity(phs)
# first_month               : first month of demand
# first_year                : first year of demand

miscellaneous = pd.read_csv("inputs_30y/miscellaneous.csv", index_col=0, header=None)
miscellaneous = miscellaneous.squeeze().copy()
"""SET HOUR BY MONTHS

Take the number of hour in the demand file.
Set up hours per months."""
first_year = miscellaneous['first_year']
first_month = miscellaneous['first_month']
first_hour = 0
last_hour = len(demand)
days_in_feb = 672
if first_year % 4 == 0:
    days_in_feb = 696

months_count = 12 * num_of_years

hours_by_months = {1: 744/time_step, 2: days_in_feb/time_step, 3: 744/time_step, 4: 720/time_step, 5: 744/time_step, 6: 720/time_step, 7: 744/time_step, 8: 744/time_step, 9: 720/time_step, 10: 744/time_step, 11: 720/time_step,
                   12: 744/time_step}

i = 1
j = first_month
hour = 0
months_hours={}
while i <= months_count:
    months_hours[i] = range(int(hour), int(hour + hours_by_months[j]))
    hour += hours_by_months[j]
    j += 1
    i += 1
    if j == 13:
        j = 1
"""SETS

Definition of set as an object of the model
"""

# Range of hour in one year
model.h = \
    pyo.RangeSet(first_hour, last_hour - 1)
# Months
model.months = \
    pyo.RangeSet(1, months_count)
# Technologies
model.tec = \
    pyo.Set(
        initialize=["onshore", "pv_g", "pv_c", "river", "lake", "biogas2", "ccgt", "nuc", "h2_ccgt",
                    "phs", "battery4", "electrolysis", "hydrogen"])
# Power plants
model.gen = \
    pyo.Set(initialize=["onshore", "pv_g", "pv_c", "river", "lake", "ccgt", "nuc"])
# Variables Technologies
model.vre = \
    pyo.Set(initialize=["onshore", "pv_g", "pv_c", "river"])
#
model.balance = \
    pyo.Set(
        initialize=["onshore", "pv_g", "pv_c", "river", "lake", "nuc", "phs", "battery4", "h2_ccgt", "ccgt"])
# Storage Technologies
model.str = \
    pyo.Set(initialize=["phs", "battery4", "hydrogen"])
# Storage Technologies
model.str_noH2 = \
    pyo.Set(initialize=["phs", "battery4"])
# Battery Storage
model.battery = \
    pyo.Set(initialize=["battery4"])

"""PARAMETERS"""

# Set the hydrogen demand for each hour
H2_demand = {}
for hour in model.h:
    H2_demand[hour] = miscellaneous['H2_demand']

"""BOUNDS VALUES

Set initial value for variables.
There is a function for each variable with bounds.
The function return the lower and the upper value.
"""


def capa_bounds(model, i):
    if i in capa_max.keys():
        return (None, capa_max[i])
    elif i == 'phs':
        return (miscellaneous['phs_discharging_lower'], miscellaneous['phs_discharging_upper'])
    else:
        return (None, None)


def s_bounds(model, i):
    if i == 'phs':
        return (miscellaneous['phs_charging_lower'], miscellaneous['phs_charging_upper'])
    else:
        return (None, None)


def capacity_bounds(model, i):
    if i == 'phs':
        return (miscellaneous['phs_energy_lower'], miscellaneous['phs_energy_upper'])
    elif i == 'hydrogen':
        return (capacity_ex['hydrogen'], None)
    else:
        return (None, None)


"""VARIABLES

Definition of variable as an object of the model
"""

# Hourly energy generation in GWh/h
model.gene = \
    pyo.Var(((tec, h) for tec in model.tec for h in model.h), within=pyo.NonNegativeReals, initialize=0)

# Hourly lost load generation in GWh/h
model.ll = \
    pyo.Var((h for h in model.h), within=pyo.NonNegativeReals, initialize=0)

# Overall yearly installed capacity in GW
model.capa = \
    pyo.Var(model.tec, within=pyo.NonNegativeReals, bounds=capa_bounds)

# Hourly electricity input of battery storage GW
model.storage = \
    pyo.Var(((storage, h) for storage in model.str for h in model.h), within=pyo.NonNegativeReals, initialize=0)

# Energy stored in each storage technology in GWh = Stage of charge
model.stored = \
    pyo.Var(((storage, h) for storage in model.str for h in model.h), within=pyo.NonNegativeReals, initialize=0)

# Charging power capacity of each storage technology
model.s = \
    pyo.Var(model.str, within=pyo.NonNegativeReals, bounds=capa_bounds)

# Energy volume of storage technology in GWh
model.capacity = \
    pyo.Var(model.str, within=pyo.NonNegativeReals, bounds=capacity_bounds)

"""FIXED VALUES"""

for tec in model.tec:
    if tec in fix_capa.keys():
        model.capa[tec].fix(fix_capa[tec])

"""CONSTRAINTS RULE

Set up a function which will return the equation of the constraint.
"""


def generation_vre_constraint_rule(model, h, vre):
    """Get constraint on variables renewable profiles generation."""

    return model.gene[vre, h] == (model.capa[vre] * load_factor[vre, h]) * time_step


def generation_capacity_constraint_rule(model, h, tec):
    """Get constraint on maximum power for non-VRE technologies."""

    return model.capa[tec] * time_step >= model.gene[tec, h]


def battery4_capacity_constraint_rule(model):
    """Get constraint on capacity of battery4."""

    return model.capa['battery4'] == model.capacity['battery4'] / 4



def combustion_2_constraint_rule(model, h):
    """Get constraint on the relationship of combustible technologies"""

    return model.gene['ccgt', h] == model.gene['biogas2', h] * miscellaneous['eta_ccgt']


def storing_constraint_rule(model, h, storage_tecs):
    """Get constraint on storing."""

    hPOne = h + 1 if h < (last_hour - 1) else 0
    charge = model.storage[storage_tecs, h] * eta_in[storage_tecs]
    discharge = model.gene[storage_tecs, h] / eta_out[storage_tecs]
    flux = charge - discharge
    return model.stored[storage_tecs, hPOne] == model.stored[storage_tecs, h] + flux


def storage_constraint_rule(model, storage_tecs):
    """Get constraint on stored energy to be equal at the end than at the start."""
    first = model.stored[storage_tecs, first_hour]
    last = model.stored[storage_tecs, last_hour - 1]
    charge = model.storage[storage_tecs, last_hour - 1] * eta_in[storage_tecs]
    discharge = model.gene[storage_tecs, last_hour - 1] / eta_out[storage_tecs]
    flux = charge - discharge
    return first == last + flux


def lake_reserve_constraint_rule(model, month):
    """Get constraint on maximum monthly lake generation."""

    return sum(model.gene['lake', hour] for hour in months_hours[month]) <= lake_inflows[month] * 1000


def stored_capacity_constraint(model, h, storage_tecs):
    """Get constraint on maximum energy that is stored in storage units"""

    return model.stored[storage_tecs, h] <= model.capacity[storage_tecs]


def storage_capacity_1_constraint_rule(model, h, storage_tecs):
    """Get constraint on the capacity with hourly charging relationship of storage"""

    return time_step * model.s[storage_tecs] >= model.storage[storage_tecs, h]


def battery_capacity_constraint_rule(model, battery):
    """Get constraint on battery's capacity."""

    return model.s[battery] == model.capa[battery]



def hydrogen_balance_constraint_rule(model, h):
    """Get constraint on hydrogen's balance."""

    gene_e_h = model.gene['electrolysis', h] + model.gene['hydrogen', h]
    dem_sto = model.gene['h2_ccgt', h] / miscellaneous['eta_h2_ccgt'] + H2_demand[h] + model.storage['hydrogen', h]
    return gene_e_h == dem_sto


def adequacy_constraint_rule(model, h):
    """Get constraint for 'supply/demand relation'"""

    sto = sum(model.storage[str_noH2, h] for str_noH2 in model.str_noH2)
    gene_electrolysis = model.gene['electrolysis', h] / miscellaneous['eta_electrolysis']
    return sum(model.gene[balance, h] for balance in model.balance) + model.ll[h] >= (
                demand[h] + sto + gene_electrolysis)


def objective_rule(model):
    """Get constraint for the final objective function."""
    fixed_cost = (sum((model.capa[tec] - capa_ex[tec]) * capex[tec] for tec in model.tec) \
                  + sum(
                (model.capacity[storage_tecs] - capacity_ex[storage_tecs]) * capex_en[storage_tecs] for storage_tecs in
                model.str) \
                  + sum(model.capa[tec] * fOM[tec] for tec in model.tec) \
                  + sum(model.s[storage_tecs] * (s_opex[storage_tecs] + s_capex[storage_tecs]) for storage_tecs in
                        model.str)) / 1000

    variable_cost = sum(sum(model.gene[tec, h] * vOM[tec] for h in model.h) for tec in model.tec) / 1000
    ll_cost = sum(model.ll[h] * 10 for h in model.h) / 1000

    return fixed_cost * num_of_years + variable_cost + ll_cost


"""CONSTRAINT CREATION

Create the constraint as an object of the model with the function declared earlier as a rule.
"""

model.generation_vre_constraint = \
    pyo.Constraint(model.h, model.vre, rule=generation_vre_constraint_rule)

model.generation_capacity_constraint = \
    pyo.Constraint(model.h, model.tec, rule=generation_capacity_constraint_rule)

model.battery_4_capacity_constraint = \
    pyo.Constraint(rule=battery4_capacity_constraint_rule)

model.combustion_2_constraint = \
    pyo.Constraint(model.h, rule=combustion_2_constraint_rule)

model.storing_constraint = \
    pyo.Constraint(model.h, model.str, rule=storing_constraint_rule)

model.storage_constraint = \
    pyo.Constraint(model.str, rule=storage_constraint_rule)

model.lake_reserve_constraint = \
    pyo.Constraint(model.months, rule=lake_reserve_constraint_rule)

model.stored_capacity_constraint = \
    pyo.Constraint(model.h, model.str, rule=stored_capacity_constraint)

model.storage_capacity_1_constraint = \
    pyo.Constraint(model.h, model.str, rule=storage_capacity_1_constraint_rule)

model.battery_capacity_constraint = \
    pyo.Constraint(model.battery, rule=battery_capacity_constraint_rule)

model.hydrogen_balance_contraint = \
    pyo.Constraint(model.h, rule=hydrogen_balance_constraint_rule)

model.adequacy_constraint = \
    pyo.Constraint(model.h, rule=adequacy_constraint_rule)

# Creation of the objective -> Cost
model.objective = pyo.Objective(rule=objective_rule)

"""SOLVE STATEMENT

Choice of the solver.
You can remove the '#' in the third line to display the output of the solver.
"""

opt = SolverFactory('gurobi')

results = opt.solve(model)

# model.display()
fixed_cost = ((sum((model.capa[tec] - capa_ex[tec]) * capex[tec] for tec in model.tec) \
               + sum(
            (model.capacity[storage_tecs] - capacity_ex[storage_tecs]) * capex_en[storage_tecs] for storage_tecs in
            model.str) \
               + sum(model.capa[tec] * fOM[tec] for tec in model.tec) \
               + sum(model.s[storage_tecs] * (s_opex[storage_tecs] + s_capex[storage_tecs]) for storage_tecs in
                     model.str)) / 1000) * num_of_years

variable_cost = sum(sum(model.gene[tec, h] * vOM[tec] for h in model.h) for tec in model.tec) / 1000

ll_cost = sum(model.ll[h] * 10 for h in model.h) / 1000
"""SET OUTPUTS VARIABLES"""

# Dictionnary which will set a little definition for each technology in the model.
technologies_definition = {
    "onshore": "onshore wind",
    "pv_g": "pv grounded",
    "pv_c": "pv commercial",
    "river": "run-of-river hydro",
    "lake": "lake and reservoirs",
    "biogas2": "biogas for ccgt",
    "ccgt": "combined cycle gas turbine",
    "nuc": "nuclear",
    "h2_ccgt": "combined cycle gas turbine using hydrogen",
    "phs": "pumped hydroelectric energy storage",
    "battery4": "4 hours battery",
    "electrolysis": "electrolysis",
    "hydrogen": "hydrogen removed from storage",
}

# The whole demand per year in TWh
sumdemand = sum(demand[hour] for hour in model.h) / 1000
# The whole electricity demand for hydrogen per year in TWh
dem_hydrogen = sum(H2_demand[hour] for hour in model.h) / 1000
# The whole generation per year in TWh
sumgene = sum(pyo.value(model.gene[gen, hour]) for hour in model.h for gen in model.gen) / 1000

# Overall yearly energy generated by the technology in TWh
gene_tec = {}
for tec in model.tec:
    gene_tec[tec] = sum(pyo.value(model.gene[tec, hour]) for hour in model.h) / 1000

    # The whole electricity input for storage per year in TWh
nSTORAGE = {}
for storage in model.str:
    for hour in model.h:
        nSTORAGE[(storage, hour)] = pyo.value(model.storage[storage, hour])

    # Electricity cost per MWh produced (euros/MWh)
lcoe_sys1 = pyo.value(model.objective) * 1000 / sumgene

# Yearly storage related loss in % of power production and in TWh
str_loss_percent = 100 * (sum(pyo.value(model.storage[storage, hour]) for storage in model.str for hour in model.h) - \
                          sum(gene_tec[storage] * 1000 for storage in model.str)) / (sumgene * 1000)
str_loss_TWh = gene_tec['electrolysis'] / miscellaneous['eta_electrolysis'] - dem_hydrogen / miscellaneous['eta_electrolysis'] - gene_tec['h2_ccgt']

for storage in model.str:
    if storage != 'hydrogen':
        str_loss_TWh += sum(nSTORAGE[storage, hour] for hour in model.h) / 1000 - gene_tec[storage]

    # Load curtailment in % of power production and in TWh
lc_percent = (100 * (sumgene - sumdemand - dem_hydrogen / 0.75) / sumgene) - str_loss_percent
lc_TWh = (sumgene - sumdemand - dem_hydrogen / 0.75) - str_loss_TWh

# Dual values
spot_price = {}
gas_price1 = {}
gas_price2 = {}
for hour in model.h:
    spot_price[hour] = - 1000000 * model.dual[model.adequacy_constraint[hour]]
    gas_price2[hour] = -1000000 * model.dual[model.combustion_2_constraint[hour]]

    # Marginal Cost
marginal_cost = sum(spot_price[hour] for hour in model.h) / (last_hour)

# Average cost of hydrogen (euros/kg)
lcoh_1 = pyo.value(model.capa['electrolysis']) * (capex['electrolysis'] + fOM['electrolysis'])
lcoh_2 = sum(pyo.value(model.gene['electrolysis', hour]) * (vOM['electrolysis'] + \
                                                            (spot_price[hour] / 1000)) for hour in model.h)
lcoh_3 = capex_en['hydrogen'] * pyo.value(model.capacity['hydrogen'])
lcoh_4 = sum(pyo.value(model.gene['electrolysis', hour]) for hour in model.h)
lcoh = (lcoh_1 + lcoh_2 + lcoh_3) * 33.33 / lcoh_4
# Electricity cost per MWh consumed (euros/MWh)
lcoe_sys2 = (pyo.value(model.objective) - (lcoh * dem_hydrogen / 33.33)) * 1000 / sumdemand

"""OUTPUTS
    There is 4 output files :
        - Summary           : A little summary with the cost and some others data
        - Hourly-Generation : Hourly data
        - Elec_Balance      : Electric Production and Consumption
        - Capacities        : List of capacities by technologies

The try, except loop is here the manage error in the creation of the outputs.
"""

# Summary

cost_dict = {}

cost_dict['total cost bEuro'] = pyo.value(model.objective) / num_of_years
cost_dict['variable cost bEuro'] = pyo.value(variable_cost) / num_of_years
cost_dict['fixed cost bEuro'] = pyo.value(fixed_cost) / num_of_years
cost_dict['ll cost'] = pyo.value(ll_cost) / num_of_years
cost_dict['lcoh'] = lcoh
cost_dict['lcoe_sys1'] = lcoe_sys1


dct1 = {k: [v] for k, v in cost_dict.items()}
df_cost = pd.DataFrame(dct1)
df_cost.to_csv(r'outputs_30y/cost_{}_{}.csv'.format(model_name, scenario))

# Hourly_Generation
df_generation = pd.DataFrame()
for tec in model.tec:
    df_generation[tec] = pyo.value(model.gene[tec, :])

df_generation['ll'] = pyo.value(model.ll[:])


for str in model.str:
    df_generation[str + 'soc'] = pyo.value(model.stored[str, :])

for str in model.str:
    df_generation[str + 'ch p '] = pyo.value(model.storage[str, :])

df_generation['demand'] =  demand


df_generation.to_csv(r"outputs_30y/generation_{}_{}.csv".format(model_name, scenario))
# price
df_price = pd.DataFrame(spot_price.values())
df_price.to_csv(r"outputs_30y/price_{}_{}.csv".format(model_name, scenario))

# Elec_balance
balance_dict = {}

for tec in model.tec:
    balance_dict[tec] = sum(pyo.value(model.gene[tec, :])) / 1000

dct1 = {k: [v] for k, v in balance_dict.items()}
df_balance = pd.DataFrame(dct1)
df_balance.to_csv(r'outputs_30y/balance_{}_{}.csv'.format(model_name, scenario))
# Capacities
capacity_dict = {}

for tec in model.tec:
    capacity_dict[tec] = pyo.value(model.capa[tec])

for str in model.str:
    capacity_dict[str + 'energy capacity'] = pyo.value(model.capacity[str])

for str in model.str:
    capacity_dict[str + 'ch p'] = pyo.value(model.s[str])

dct_capa = {k: [v] for k, v in capacity_dict.items()}  # WORKAROUND
df_capa = pd.DataFrame(dct_capa)
df_capa.to_csv(r'outputs_30y/capacity_{}_{}.csv'.format(model_name, scenario))