# An economic model/simulation to back up my arguments here
# http://vitalik.ca/jekyll/general/2017/03/11/a_note_on_charity.html

# Marginal cost of creating a sandwich
sandwich_production_cost = 6

# A list of agents, each of the form (sandwich_reserve_price, charity_utility),
# where sandwich_reserve_price is the net utility in dollars they get from a
# sandwich and charity_utility is the net utility in dollars (including
# feelings of warm altruistic fuzziness) they get from volunteering at a
# charity.
people = []
for i in range(301):
    for j in range(301):
        # Sandwich reserve price evenly distributed between $6 and $9
        # Charity utility evenly distributed between -$4.5 and $1.5
        people.append((6.0001 + i * 0.01, -4.501 + j * 0.02))

people_sorted_by_charity_utility = sorted(people, key=lambda f: f[1])

# Get the profit a merchant gets and the number of charity volunteers, given a
# particular price that the merchant assigns to charity volunteers and to
# non-volunteers
def get_profit_and_charity_goers(sandwich_no_charity_price, sandwich_charity_price):

    charity_goers = 0
    sandwich_buyers = 0
    revenue = 0
    
    for (sandwich_reserve_price, charity_utility) in people:
        # Utility to this agent from doing nothing
        u_do_nothing = 0
        # Utility from buying a sandwich
        u_buy_sandwich = sandwich_reserve_price - sandwich_no_charity_price
        # Utility from going to charity
        u_charity = charity_utility
        # Utility from buying a sandwich and going to charity
        u_charity_plus_sandwich = charity_utility + sandwich_reserve_price - sandwich_charity_price
        # Max utility (if any utility equals the max utility that means
        # that action is taken)
        u_max = max(u_do_nothing, u_buy_sandwich, u_charity, u_charity_plus_sandwich)
        # Does this agent volunteer at the charity?
        if u_max in (u_charity, u_charity_plus_sandwich):
            charity_goers += 1
        # Does this agent buy a sandwich?
        if u_max in (u_buy_sandwich, u_charity_plus_sandwich):
            sandwich_buyers += 1
            revenue += sandwich_charity_price if u_charity_plus_sandwich == u_max else sandwich_no_charity_price
    
    # Profit = revenue minus cost of production
    profit = revenue - sandwich_buyers * sandwich_production_cost
    return profit, charity_goers

# Compare non-price-discriminating case with various degrees of
# discounts for charity volunteers
p1, c1 = get_profit_and_charity_goers(7.5, 7.5)
for i in range(6):
    p2, c2 = get_profit_and_charity_goers(7.5 + i * 0.02, 7.5 - i * 0.06)
    print('At prices (%.2f, %.2f)' % (7.5 + i * 0.02, 7.5 - i * 0.06))
    print('Net cost to merchant: %.2f' % (p1 - p2))
    print('Net gain in charity goers: %.2f' % (c2 - c1))
    # Calculate the subsidy that would lead to an equal increase in charity
    # contributors
    equiv_subsidy_size = people_sorted_by_charity_utility[-c1][1] - \
                         people_sorted_by_charity_utility[-c2][1]
    # The subsidy would go to charity workers only
    print('Cost of providing equipotent straight subsidy: %.2f'
           % (equiv_subsidy_size * c2))

#c2 * i * 0.08)
