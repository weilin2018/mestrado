

class animation(object):

    def __init__(self,kind='elev',shape='2D'):
        if shape=='3D':
            # for 3D-parameters, such as temperature, salinity, velocity
            self.var = self.ncin[kind][self.timeStart:self.timeEnd,:,:,:].values
        elif shape=='2D':
            # elevation
            self.var = self.ncin[kind][self.timeStart:self.timeEnd,:,:].values
