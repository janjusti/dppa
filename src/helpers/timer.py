import time 

class CustomTimer:
    def __init__(self):
        self.init_time = time.time()
        self.prev_time = self.init_time
        print('+++ Timer started.')
    
    def set_cp(self):
        self.curr_time = time.time()
        cp_time = self.curr_time - self.prev_time
        self.prev_time = self.curr_time        
        print('--- Lap time:', str(round(cp_time, 3)) + 's')
    
    def end(self):
        total_time = time.time() - self.init_time
        print("<<<<<<<< END >>>>>>>>")
        print('+++ Execution time:', str(round(total_time, 3)) + 's')
        print("<<<<<<<< END >>>>>>>>")