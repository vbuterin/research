# An economic model/simulation to back up my arguments here
# http://vitalik.ca/jekyll/general/2017/03/11/a_note_on_charity.html
 
# Marginal cost of creating a sandwich
sandwich_production_cost = 6

PRECISION = 10
 
# A list of agents, each of the form (sandwich_reserve_price, charity_utility),
# where sandwich_reserve_price is the net utility in dollars they get from a
# sandwich and charity_utility is the net utility in dollars (including
# feelings of warm altruistic fuzziness) they get from volunteering at a
# charity.
people = []
# Total number of people is 9025, which is close to the 9000 from Vitalik's blog post
for i in range(95 * PRECISION):
    for j in range(95 * PRECISION):
        # Sandwich reserve price evenly distributed between $6 and $9
        # Charity utility evenly distributed between -$2.25 and $.75
        # This gives both distributions the same slope to make the
        # numbers nicely line up with Vitalik's O(d^2) arguments
        people.append((6.00 + i * 3./(95 * PRECISION - 1), -2.25 + j * 3./(95 * PRECISION - 1)))
 
people_sorted_by_charity_utility = sorted(people, key=lambda f: f[1])
 
# Get the profit a merchant gets and the number of charity volunteers, given a
# particular price that the merchant assigns to charity volunteers and to
# non-volunteers
def get_demand_profit_and_charity_goers(sandwich_no_charity_price, sandwich_charity_price):
 
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
    return sandwich_buyers / PRECISION**2, profit / PRECISION**2, charity_goers / PRECISION**2
 
# Compare non-price-discriminating case with various degrees of
# discounts for charity volunteers
d1, p1, c1 = get_demand_profit_and_charity_goers(7.5, 7.5)
for (hp1, hp2) in [(7.5 + i * .05, 7.5 - i * .15) for i in range(1,2)]:
  d2, p2, c2 = get_demand_profit_and_charity_goers(hp1, hp2)
  print('Original demand: %d' % d1)
  print('New demand: %d' % d2)
  print('Original profit: %.2f' % p1)
  print('New profit: %.2f' % p2)
  print('At prices (%.2f, %.2f)' % (hp1, hp2))
  print('Net cost to merchant: $%.2f (much higher than $22.50, edit blog?)' % (p1 - p2))
  print('Net gain in charity goers: %d new workers' % (c2 - c1))
  # Calculate the subsidy that would lead to an equal increase in charity
  # contributors
  equiv_subsidy_size = people_sorted_by_charity_utility[-c1*PRECISION**2][1] - \
                       people_sorted_by_charity_utility[-c2*PRECISION**2][1]
  # The subsidy would go to charity workers only
  print('Cost of providing equipotent straight subsidy: $%.2f'
         % (equiv_subsidy_size * c2))
  # Cost to charity of providing all charity workers with $.20 hoagie subsidies
  # This is more efficient than direct subsidy for the same reason that cable
  # providers offer package deals with a bunch of stuff nobody actually wants
  print('Cost of providing equipotent all-worker hoagie promotion: $%.2f'
         % (.2 * c2 * (9-hp2)/3))
  # Cost to charity of providing new-hires with $.20 hoagie subsidies
  # Optimal price discrimination between new hires and legacy workers
  # An even better strategy would be to pay each worker exactly their reserve price
  # but they weren't doing that already, so we'll just go with this
  print('Cost of providing equipotent new-hire hoagie promotion: $%.2f'
         % (.2 * (c2-c1)))
  # Regular donations are more versatile and therefore preferred
  print('Lower bound on deadweight loss of sticker strategy: $%.2f'
         % ((p1-p2) - .2 * (c2-c1)))
