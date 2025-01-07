from kaggle_environments import make
env = make("chess", debug=True)

# Define the number of iterations to run
num_iterations = 2

# Initialize variables to store the cumulative results
mainscore = 0
main2score = 0


for i in range(num_iterations):
    l = 0
    # Run the first match: main.py vs main2.py
    result = env.run(["main.py", "main2.py"])
    print(f"Iteration {i+1}: main.py vs main2.py")
    print("Agent exit status/reward/time left: ")
    for agent in result[-1]:
        l += 1
        print("\t", agent.status, "/", agent.reward, "/", agent.observation.remainingOverageTime)
        if l == 1:
            mainscore += agent.reward
        elif l == 2:
            main2score += agent.reward
        else:
            print("error no marks1")
    print("\n")
    l = 0
    env.render(mode="ipython", width=600, height=600)

    result = env.run(["main2.py", "main.py"])
    print(f"Iteration {i+1}: main2.py vs main.py")
    print("Agent exit status/reward/time left: ")
    for agent in result[-1]:
        l += 1
        print("\t", agent.status, "/", agent.reward, "/", agent.observation.remainingOverageTime)
        if l == 1:
            main2score += agent.reward
        elif l == 2:
            mainscore += agent.reward
        else:
            print("error no marks2")
    print("\n")
    
    env.render(mode="ipython", width=600, height=600)
    print("SCORES: ", "main.py: ",mainscore, " main2.py: ",main2score)




