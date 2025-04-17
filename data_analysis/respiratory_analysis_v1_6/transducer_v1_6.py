class Transducer:

    def __init__(self, supply_voltage, min_inches_water, max_inches_water, ref_pascals):
        """Initialize the PressureTransducer object
           supply_voltage: voltage supplied to the transducer in volts
           min_inches_water: transducer's minimum pressure rating in inches of H2O
           max_inches_water: transducer's maximum pressure rating in inches of H2O
           ref_pascals: pressure supplied to the transducer's reference port in Pascals"""
        self.supply_voltage = supply_voltage
        self.min_inches_water = min_inches_water
        self.max_inches_water = max_inches_water
        self.ref_pascals = ref_pascals
        self.min_pascals = inches_water_to_pascals(min_inches_water, ref_pascals)
        self.max_pascals = inches_water_to_pascals(max_inches_water, ref_pascals)


def inches_water_to_pascals(inches_water, ref_inches_water):
    conversion_factor = 248.8  # Pascals / inch H2O
    pascals = ref_inches_water + (inches_water * conversion_factor)
    return pascals


def main():
    supply_voltage = 5.0
    min_inches_water = -4.0
    max_inches_water = 4.0
    ref_pascals = 93073  # Atmospheric pressure [Pascals] for Tucson, AZ at 22C and 728m altitude
    my_transducer = Transducer(supply_voltage, min_inches_water, max_inches_water, ref_pascals)
    print(my_transducer.min_pascals, my_transducer.max_pascals)


if __name__ == "__main__":
    main()
