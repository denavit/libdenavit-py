
class AnalysisResults:
    print_each_analysis_time_increment = True
    total_analysis_time = 0.
    
    def __init__(self):
        pass
        
    def add_to_analysis_time(self,tic,toc):
        self.total_analysis_time += toc - tic
        if self.print_each_analysis_time_increment:
            print(f"Adding {toc - tic:0.4f} seconds to total analysis run time.")

    def print_total_analysis_time(self):
        print(f"Total analysis time: {self.total_analysis_time:0.4f} seconds")
