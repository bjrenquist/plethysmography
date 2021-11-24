import respiration_v1_5 as respiration
import transducer_v1_5 as transducer
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.widgets import Slider


class AdjustAxisDialog(simpledialog.Dialog):
    def __init__(self, parent, axes, x_or_y, axis_min, axis_max, title=None):
        self.axes = axes
        self.x_or_y = x_or_y
        self.axis_min = axis_min
        self.axis_max = axis_max
        self.entry_min = None
        self.entry_max = None
        self.entry_min_value = None
        self.entry_max_value = None
        self.label_min = None
        self.label_max = None
        simpledialog.Dialog.__init__(self, parent, title)

    def body(self, parent):
        self.entry_min_value = tk.StringVar()
        self.entry_max_value = tk.StringVar()
        self.entry_min = tk.Entry(parent, textvariable=self.entry_min_value)
        self.entry_min_value.set(self.axis_min)
        self.entry_max = tk.Entry(parent, textvariable=self.entry_max_value)
        self.entry_max_value.set(self.axis_max)
        self.label_min = tk.Label(parent, text='Minimum Value:')
        self.label_max = tk.Label(parent, text='Maximum Value:')
        self.label_max.grid(sticky=tk.W, padx=5)
        self.entry_max.grid(sticky=tk.W + tk.E, padx=5)
        self.label_min.grid(sticky=tk.W, padx=5)
        self.entry_min.grid(sticky=tk.W + tk.E, padx=5)

    def apply(self):
        axis_min = float(self.entry_min.get())
        axis_max = float(self.entry_max.get())
        for idx, axis in enumerate(self.axes):
            if self.x_or_y[idx] == 'x':
                axis.set_xlim(axis_min, axis_max)
                self.parent.fig.canvas.draw_idle()
            elif self.x_or_y[idx] == 'y':
                axis.set_ylim(axis_min, axis_max)
                self.parent.fig.canvas.draw_idle()
            else:
                # print('choose x or y axis to adjust')
                return


class AnalysisWindow(tk.Frame):

    def __init__(self, parent=None):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.grid()

        self.title = 'Respiration Analysis'
        self.winfo_toplevel().title(self.title)

        # Todo: update self.version for new software versions
        self.version = 'v1_5'
        # version 1_5 updates the autocalibration step width threshold parameter to default to 0.5s - the previous
        # default of 2.0s was too long to pick up shorter calibration injections, and the decrease does not appear to
        # effect calibration negatively. The global drift adjustment was reworked to avoid strange linear algebra
        # library error involving scipy.detrend(). Some improvements were made to the legend of the calibration
        # injection points regression plot; the number of included and excluded points are now reported.

        self.filename = ''
        self.short_filename = ''
        self.use_prev_cal = False
        self.prev_cal_filename = ''
        self.data = None

        self.widget_padding = 10
        self.start_time = 0.0
        self.stop_time = 20.0
        self.time_window = 10.0
        self.roi_start = self.start_time
        self.roi_start_idx = 0
        self.roi_stop = self.stop_time
        self.roi_stop_idx = 1

        # Todo: Allow user to select different transducers
        # Transducer settings
        self.supply_voltage = 5.0
        self.min_voltage = 0.0
        self.max_voltage = 5.0
        self.min_inches_water = -4.0
        self.max_inches_water = 4.0
        self.ref_pressure = 93073  # Atmospheric pressure in Pascals for Tucson, AZ at 22C and 728m altitude
        self.transducer = transducer.Transducer(self.supply_voltage, self.min_inches_water, self.max_inches_water,
                                                self.ref_pressure)
        self.min_pressure = self.transducer.min_pascals  # min transducer pressure in Pascals
        self.max_pressure = self.transducer.max_pascals  # max transducer pressure in Pascals

        self.sat_min_pressure = float()
        self.sat_max_pressure = float()
        self.axis_min_pressure = float()
        self.axis_max_pressure = float()

        self.min_injection_volume = 0.0
        self.max_injection_volume = 0.5
        self.min_volume = -0.5
        self.max_volume = 0.5
        self.sat_min_volume = float()
        self.sat_max_volume = float()
        self.min_flow = -20.0
        self.max_flow = 20.0

        # Breathing analysis algorithm settings
        self.sat_min_voltage = 0.150  # lower saturation voltage
        self.sat_max_voltage = 4.800  # lower saturation voltage
        self.volume_calibration_cutoff = 0.5  # previously 0.5
        self.volume_calibration_order = 6
        self.volume_smoothing_sigma = 5.0
        self.autocal_sigma = 0.4
        self.autocal_search_left = 24
        self.autocal_search_right = 12
        self.autocal_dup_comment_distance = 2
        self.autocal_lowest_cal_vol = 0.05
        self.autocal_highest_cal_vol = 0.5
        self.autocal_ref_cal_vol = 0.1
        self.autocal_w_ref_cal_vol = 2
        self.autocal_w_variance = 0.4
        self.autocal_num_std_step_height_threshold = 0.25
        self.autocal_step_width_threshold = 0.5
        self.autocal_lower_variance_ratio = 0.5
        self.autocal_downsampling_rate = 4
        self.autocal_num_std_exclude = 1.0
        self.autocal_make_plots = False

        self.local_drift_correction_window_size = 2.0
        self.flow_oversmoothing_sigma = 30.0
        self.peak_detection_window_sizes = [10, 20, 50, 100, 200]
        self.peak_detection_shift_sizes = [0.000, 0.333, 0.667]
        self.pause_search_threshold = 25.0  # Initial percent of inhale/exhale volume during which to look for a pause
        self.number_histogram_bins = 50
        self.extra_zero_bin_threshold = 0.10
        self.extra_pause_bin_threshold = 0.10
        self.mode_bin_ratio = 2
        self.extra_pause_bin_ratio = 0.75

        self.movement_threshold = 0.4

        self.exclusion_adjacency = 3
        self.must_be_between_exclusions = True
        self.avg_window_size = 60.0
        self.avg_window_offset = 30.0

        self.outlier_test_stats = ['piv',
                                   'pev',
                                   'pif',
                                   'pef',
                                   'ti',
                                   'te']
        self.outlier_window_size = 30.0  # seconds, centered at current breath's start_idx
        self.max_num_stdevs = 4.0
        self.graphical_downsampling_rate = 1
        self.display_smoothed_flow = tk.BooleanVar()
        self.display_smoothed_flow.set(False)
        self.save_signal_data = tk.BooleanVar()
        self.save_signal_data.set(False)
        self.remove_saturated = tk.BooleanVar()
        self.remove_saturated.set(True)
        self.remove_movement = tk.BooleanVar()
        self.remove_movement.set(True)
        self.remove_outlier = tk.BooleanVar()
        self.remove_outlier.set(True)

        # ax1 for voltage over time
        self.ax1_xmin = self.start_time
        self.ax1_xmax = self.start_time + self.time_window
        self.ax1_ymin = self.min_voltage
        self.ax1_ymax = self.max_voltage
        self.ax1_xlabel = 'Time [s]'
        self.ax1_ylabel = 'Voltage [V]'

        # ax2 for pressure over time
        self.ax2_xmin = self.start_time
        self.ax2_xmax = self.start_time + self.time_window
        self.ax2_ymin = self.min_pressure
        self.ax2_ymax = self.max_pressure
        self.ax2_xlabel = 'Time [s]'
        self.ax2_ylabel = 'Pressure [Pa]'

        # ax3 for pressure over injection volume
        self.ax3_xmin = self.min_injection_volume
        self.ax3_xmax = self.max_injection_volume
        self.ax3_ymin = 0
        self.ax3_ymax = self.max_pressure - self.ref_pressure
        self.ax3_xlabel = 'Change in Chamber Volume [mL]'
        self.ax3_ylabel = 'Change in Pressure [Pa]'

        # ax4 for volume over time
        self.ax4_xmin = self.start_time
        self.ax4_xmax = self.start_time + self.time_window
        self.ax4_ymin = self.min_volume
        self.ax4_ymax = self.max_volume
        self.ax4_xlabel = 'Time [s]'
        self.ax4_ylabel = 'Volume [mL]'

        # ax5 for flow over time
        self.ax5_xmin = self.start_time
        self.ax5_xmax = self.start_time + self.time_window
        self.ax5_ymin = self.min_flow
        self.ax5_ymax = self.max_flow
        self.ax5_xlabel = 'Time [s]'
        self.ax5_ylabel = 'Flow [mL/s]'

        # ax6 for overlaid breaths
        self.ax6_xlabel = 'Time [s]'
        self.ax6_ylabel = 'Volume [mL]'

        # Configure menu bar
        self.label_file_menu = 'File'
        self.label_open = 'Open...'
        self.label_save = 'Save As...'
        self.label_exit = 'Exit'

        self.label_edit_menu = 'Edit'
        self.label_sat_min_voltage = 'Set Lower Saturation Voltage'
        self.label_sat_max_voltage = 'Set Upper Saturation Voltage'
        self.label_filter_cutoff = 'Set Lowpass Filter Cutoff Frequency'
        self.label_filter_order = 'Set Lowpass Filter Order'
        self.label_autocal_settings = 'Autocalibration Settings'
        self.label_volume_smoothing_sigma = 'Set Volume Signal Gaussian Smoothing Filter Sigma'
        self.label_local_drift_corr_window_size = 'Set Local Drift Correction Window Size'
        self.label_flow_oversmoothing_sigma = 'Set Flow Signal Gaussian Smoothing Filter Sigma'
        self.label_peak_detection_window_sizes = 'Set Peak Detection Window Sizes'
        self.label_peak_detection_shift_sizes = 'Set Peak Detection Shift Sizes'
        self.label_pause_search_threshold = 'Set Pause Search Volume Threshold'
        self.label_number_histogram_bins = 'Set Pause Search Number of Histogram Bins'
        self.label_extra_zero_bin_threshold = 'Set Pause Search Extra Zero Bin Threshold'
        self.label_extra_pause_bin_threshold = 'Set Pause Search Extra Pause Bin Threshold'
        self.label_mode_bin_ratio = 'Set Pause Search Mode Bin Ratio'
        self.label_extra_pause_bin_ratio = 'Set Pause Search Extra Pause Bin Ratio'
        self.label_movement_threshold = 'Set Movement Threshold'
        self.label_exclusion_adjacency = 'Set Breath Exclusion Adjacency'
        self.label_must_be_between_exlusions = 'Set Double-Sided Exclusion Adjacency Requirement'
        self.label_avg_window_size = 'Set Averaging Window Size'
        self.label_avg_window_offset = 'Set Averaging Window Offset'
        self.label_graphical_downsampling_rate = 'Set Graphical Downsampling Rate'
        self.label_outlier_settings = 'Outlier Settings'
        self.label_display_smoothed_flow = 'Display Over-Smoothed Flow Signal?'
        self.label_save_signal_data = 'Save Signal Data?'
        self.label_remove_saturated = 'Remove Breaths with Voltage Saturation?'
        self.label_remove_movement = 'Remove Breaths where Movement Detected?'
        self.label_remove_outlier = 'Remove Outlier Breaths?'

        self.label_analysis_menu = 'Analysis'
        self.label_start_auto_calibration = 'Start Automatic Volume Calibration'
        self.label_start_volume_calibration = 'Start Manual Volume Calibration'
        self.label_stop_volume_calibration = 'Stop Manual Volume Calibration'
        self.label_find_breath_landmarks = 'Find Breath Landmarks'

        self.label_view_menu = 'View'
        self.label_time_scale = 'Adjust Time Axis Scale'
        self.label_go_to_time = 'Go To Time...'
        self.label_voltage_scale = 'Adjust Voltage Axis Scale'
        self.label_pressure_scale = 'Adjust Pressure Axis Scale'
        self.label_inj_volume_scale = 'Adjust Calibration Change in Chamber Volume Scale'
        self.label_cal_pressure_scale = 'Adjust Calibration Change in Pressure Scale'
        self.label_volume_scale = 'Adjust Volume Axis Scale'
        self.label_flow_scale = 'Adjust Flow Axis Scale'

        self.label_help_menu = 'Help'

        # menus dict - 'Menu Title': tuple(['Command Title', command, state])
        self.menus = {self.label_file_menu: (
                      ([tk.COMMAND, self.label_open, self.open, tk.NORMAL]),
                      ([tk.COMMAND, self.label_save, self.save, tk.DISABLED]),
                      ([tk.SEPARATOR]),
                      ([tk.COMMAND, self.label_exit, self.close_window, tk.NORMAL])),

                      self.label_edit_menu: (
                      ([tk.COMMAND, self.label_sat_min_voltage, self.set_sat_min_voltage, tk.NORMAL]),
                      ([tk.COMMAND, self.label_sat_max_voltage, self.set_sat_max_voltage, tk.NORMAL]),
                      ([tk.COMMAND, self.label_filter_cutoff, self.set_filter_cutoff, tk.NORMAL]),
                      ([tk.COMMAND, self.label_filter_order, self.set_filter_order, tk.NORMAL]),
                      ([tk.COMMAND, self.label_autocal_settings, self.set_autocal_settings, tk.NORMAL]),
                      ([tk.COMMAND, self.label_volume_smoothing_sigma, self.set_smoothing_sigma, tk.NORMAL]),
                      ([tk.COMMAND, self.label_local_drift_corr_window_size, self.set_local_drift_corr_window_size,
                        tk.NORMAL]),
                      ([tk.COMMAND, self.label_flow_oversmoothing_sigma, self.set_flow_oversmoothing_sigma, tk.NORMAL]),
                      ([tk.COMMAND, self.label_peak_detection_window_sizes, self.set_peak_detection_window_sizes,
                        tk.NORMAL]),
                      ([tk.COMMAND, self.label_peak_detection_shift_sizes, self.set_peak_detection_shift_sizes,
                        tk.NORMAL]),
                      ([tk.COMMAND, self.label_pause_search_threshold, self.set_pause_search_threshold, tk.NORMAL]),
                      ([tk.COMMAND, self.label_number_histogram_bins, self.set_number_histogram_bins, tk.NORMAL]),
                      ([tk.COMMAND, self.label_extra_zero_bin_threshold, self.set_extra_zero_bin_threshold, tk.NORMAL]),
                      ([tk.COMMAND, self.label_extra_pause_bin_threshold, self.set_extra_pause_bin_threshold,
                        tk.NORMAL]),
                      ([tk.COMMAND, self.label_mode_bin_ratio, self.set_mode_bin_ratio, tk.NORMAL]),
                      ([tk.COMMAND, self.label_extra_pause_bin_ratio, self.set_extra_pause_bin_ratio, tk.NORMAL]),
                      ([tk.COMMAND, self.label_movement_threshold, self.set_movement_threshold, tk.NORMAL]),
                      ([tk.COMMAND, self.label_exclusion_adjacency, self.set_exclusion_adjacency, tk.NORMAL]),
                      ([tk.COMMAND, self.label_must_be_between_exlusions, self.set_must_be_between_exclusions,
                        tk.NORMAL]),
                      ([tk.COMMAND, self.label_avg_window_size, self.set_avg_window_size, tk.NORMAL]),
                      ([tk.COMMAND, self.label_avg_window_offset, self.set_avg_window_offset, tk.NORMAL]),
                      ([tk.COMMAND, self.label_outlier_settings, self.set_outlier_settings, tk.NORMAL]),
                      ([tk.COMMAND, self.label_graphical_downsampling_rate, self.set_downsampling_rate, tk.NORMAL]),
                      ([tk.SEPARATOR]),
                      ([tk.CHECKBUTTON, self.label_display_smoothed_flow, True, False, self.display_smoothed_flow]),
                      ([tk.CHECKBUTTON, self.label_remove_outlier, True, False, self.remove_outlier]),
                      ([tk.CHECKBUTTON, self.label_remove_saturated, True, False, self.remove_saturated]),
                      ([tk.CHECKBUTTON, self.label_remove_movement, True, False, self.remove_movement]),
                      ([tk.CHECKBUTTON, self.label_save_signal_data, True, False, self.save_signal_data])),

                      self.label_analysis_menu: (
                      ([tk.COMMAND, self.label_start_auto_calibration, self.start_auto_calibration, tk.DISABLED]),
                      ([tk.SEPARATOR]),
                      ([tk.COMMAND, self.label_start_volume_calibration, self.start_volume_calibration, tk.DISABLED]),
                      ([tk.COMMAND, self.label_stop_volume_calibration, self.stop_volume_calibration, tk.DISABLED]),
                      ([tk.SEPARATOR]),
                      ([tk.COMMAND, self.label_find_breath_landmarks, self.find_breath_landmarks, tk.DISABLED])),

                      self.label_view_menu: (
                      ([tk.COMMAND, self.label_go_to_time, self.go_to_time, tk.DISABLED]),
                      ([tk.COMMAND, self.label_time_scale, self.adjust_time_scale, tk.DISABLED]),
                      ([tk.COMMAND, self.label_voltage_scale,
                        lambda: AdjustAxisDialog(self, [self.ax1], ['y'], self.ax1_ymin, self.ax1_ymax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]),
                      ([tk.COMMAND, self.label_pressure_scale,
                        lambda: AdjustAxisDialog(self, [self.ax2], ['y'], self.ax2_ymin, self.ax2_ymax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]),
                      ([tk.COMMAND, self.label_inj_volume_scale,
                        lambda: AdjustAxisDialog(self, [self.ax3], ['x'], self.ax3_xmin, self.ax3_xmax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]),
                      ([tk.COMMAND, self.label_cal_pressure_scale,
                        lambda: AdjustAxisDialog(self, [self.ax3], ['y'], self.ax3_ymin, self.ax3_ymax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]),
                      ([tk.COMMAND, self.label_volume_scale,
                        lambda: AdjustAxisDialog(self, [self.ax4], ['y'], self.ax4_ymin, self.ax4_ymax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]),
                      ([tk.COMMAND, self.label_flow_scale,
                        lambda: AdjustAxisDialog(self, [self.ax5], ['y'], self.ax5_ymin, self.ax5_ymax,
                                                 title='Voltage Axis Settings'), tk.DISABLED]))
                      }

        # TODO: Add help options; open readme files?
        # 'Help': tuple([1, 3, 4])}

        self.menu_bar = tk.Menu(self.parent)
        self.submenus = {}
        for key, values in self.menus.items():
            menu = tk.Menu(self.menu_bar, tearoff=0)
            for menu_type, *params in values:
                if menu_type == tk.SEPARATOR:
                    menu.add_separator()
                elif menu_type == tk.COMMAND:
                    label, command, state = params
                    menu.add_command(label=label, command=command, state=state)
                elif menu_type == tk.CHECKBUTTON:
                    label, onvalue, offvalue, variable = params
                    menu.add_checkbutton(label=label, onvalue=onvalue, offvalue=offvalue, variable=variable)
            self.submenus[key] = menu
            self.menu_bar.add_cascade(label=key, menu=menu)
        self.parent.config(menu=self.menu_bar)

        # Configure plot area
        self.fig = plt.Figure()
        # self.fig2 = plt.figure()
        self.fig.set_size_inches(10.0, 6.0)
        self.fig.subplots_adjust(bottom=0.05, top=0.95, hspace=0.5)
        self.spec = gridspec.GridSpec(ncols=1, nrows=3, height_ratios=[10, 10, 1], figure=self.fig)
        self.canvas = FigureCanvasTkAgg(self.fig, self.parent)
        self.canvas.get_tk_widget().grid()

        self.ax1 = self.fig.add_subplot(self.spec[0, 0], label='1')
        self.ax1.set_xlim(self.ax1_xmin, self.ax1_xmax)
        self.ax1.set_ylim(self.ax1_ymin, self.ax1_ymax)
        self.ax1.set_xlabel(self.ax1_xlabel)
        self.ax1.set_ylabel(self.ax1_ylabel)
        self.ax2 = self.fig.add_subplot(self.spec[1, 0], label='2')
        self.ax2.set_xlim(self.ax2_xmin, self.ax2_xmax)
        self.ax2.set_ylim(self.ax2_ymin, self.ax2_ymax)
        self.ax2.set_xlabel(self.ax2_xlabel)
        self.ax2.set_ylabel(self.ax2_ylabel)
        self.plot_selection = None
        self.ax3 = self.fig.add_subplot(self.spec[0, 0], label='3')
        self.ax3_line1, = self.ax3.plot([], [], 'ko', picker=5)
        self.fig.canvas.mpl_connect('pick_event', self.edit_calibration_point)
        self.ax3_line2, = self.ax3.plot([], [], 'k:')
        self.ax3.set_xlim(self.ax3_xmin, self.ax3_xmax)
        self.ax3.set_ylim(self.ax3_ymin, self.ax3_ymax)
        self.ax3.set_xlabel(self.ax3_xlabel)
        self.ax3.set_ylabel(self.ax3_ylabel)
        self.ax4 = self.fig.add_subplot(self.spec[0, 0], label='4')
        self.ax4.set_xlim(self.ax4_xmin, self.ax4_xmax)
        self.ax4.set_ylim(self.ax4_ymin, self.ax4_ymax)
        self.ax4.set_xlabel(self.ax4_xlabel)
        self.ax4.set_ylabel(self.ax4_ylabel)
        self.ax5 = self.fig.add_subplot(self.spec[1, 0], label='5')
        self.ax5.set_xlim(self.ax5_xmin, self.ax5_xmax)
        self.ax5.set_ylim(self.ax5_ymin, self.ax5_ymax)
        self.ax5.set_xlabel(self.ax5_xlabel)
        self.ax5.set_ylabel(self.ax5_ylabel)
        # self.ax6 = self.fig2.add_subplot(1, 1, 1)
        # self.ax6.set_xlabel(self.ax6_xlabel)
        # self.ax6.set_ylabel(self.ax6_ylabel)
        for ax in [self.ax3, self.ax4, self.ax5]:
            self.fig.delaxes(ax)

        self.axpos = self.fig.add_subplot(self.spec[2, 0], label='xpos')
        self.axpos.set_xlabel('\u2190 X Position \u2192')
        self.spos = Slider(self.axpos, '', self.start_time,
                           self.stop_time - self.time_window,
                           valinit=self.start_time)

        def update(pos):
            for axes in [self.ax1, self.ax2, self.ax4, self.ax5]:
                axes.set_xlim(pos, pos + self.time_window)
            self.fig.canvas.draw_idle()

        self.spos.on_changed(update)

    def adjust_time_scale(self):
        new_time_window = simpledialog.askfloat(self.label_time_scale, 'Seconds to show')
        if not new_time_window:
            return
        self.time_window = new_time_window
        for ax in [self.ax1, self.ax2, self.ax4, self.ax5]:
            ax.set_xlim(self.spos.val, self.spos.val + self.time_window)

        self.spos.valmin = self.start_time
        self.spos.valmax = self.stop_time - self.time_window
        self.spos.ax.set_xlim(self.spos.valmin, self.spos.valmax)

        self.ax1_xmin, self.ax1_xmax = self.ax1.get_xlim()
        self.ax2_xmin, self.ax2_xmax = self.ax2.get_xlim()
        self.ax4_xmin, self.ax1_xmax = self.ax4.get_xlim()
        self.ax5_xmin, self.ax2_xmax = self.ax5.get_xlim()

        self.fig.canvas.draw_idle()

    def calculate_volumetric_flow(self):
        # Calculate and plot flow signal
        self.data.calculate_flow_signal(self.data.local_corrected_volume)
        self.hide_axes(self.fig, self.ax2)
        self.show_axes(self.fig, self.ax5)
        self.plot_signal(self.ax5, self.data.time, self.data.flow, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=True, color='k', linestyle='-')

        # Update menu bar options
        self.update_menu_bar_state(self.label_analysis_menu, self.label_find_breath_landmarks, tk.NORMAL)
        self.update_menu_bar_state(self.label_view_menu, self.label_pressure_scale, tk.DISABLED)
        self.update_menu_bar_state(self.label_view_menu, self.label_flow_scale, tk.NORMAL)

    def check_calibration_axes(self):
        # Axis min/max volume and pressure values
        volume_min = self.ax3_xmin
        volume_max = self.ax3_xmax
        pressure_min = self.ax3_ymin
        pressure_max = self.ax3_ymax

        # Find min/max volume and pressure values from all calibration points
        cal_volumes = self.data.calibration_volumes
        cal_pressures = self.data.calibration_pressures

        cal_volumes_min = np.amin(cal_volumes)
        cal_volumes_max = np.amax(cal_volumes)
        cal_pressures_min = np.amin(cal_pressures)
        cal_pressures_max = np.amax(cal_pressures)

        # Check if calibration point min/max values lie within min/max boundaries of axes
        scaling_factor = 0.90
        if cal_volumes_min < scaling_factor * volume_min:
            # Adjust volume axis min boundary
            new_volume_min = np.floor(10.0 * cal_volumes_min / scaling_factor) / 10.0
            self.ax3_xmin = new_volume_min

        if cal_volumes_max > scaling_factor * volume_max:
            # Adjust volume axis max boundary
            new_volume_max = np.ceil(10.0 * cal_volumes_max / scaling_factor) / 10.0
            self.ax3_xmax = new_volume_max

        if cal_pressures_min < scaling_factor * pressure_min:
            # Adjust pressure axis min boundary
            new_pressure_min = np.floor(10.0 * cal_pressures_min / scaling_factor) / 10.0
            self.ax3_ymin = new_pressure_min

        if cal_pressures_max > scaling_factor * pressure_max:
            # Adjust pressure axis max boundary
            new_pressure_max = np.ceil(10.0 * cal_pressures_max / scaling_factor) / 10.0
            self.ax3_ymax = new_pressure_max

        # Redraw axes
        self.ax3.set_xlim(self.ax3_xmin, self.ax3_xmax)
        self.ax3.set_ylim(self.ax3_ymin, self.ax3_ymax)
        self.fig.canvas.draw_idle()

        return

    def close_window(self):
        # Properly close out matplotlib figures and destroy tkinter root window on closing window with 'x'
        if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            plt.close('all')
            self.master.destroy()

    def display_breaths(self):
        # Todo: Allow user to select which landmarks to display in settings window
        # Setup arrays for volume and flow extrema
        inhale_flow_peak_t = np.full(self.data.num_breaths, np.nan)
        inhale_flow_peak_f = np.full(self.data.num_breaths, np.nan)

        exhale_flow_peak_t = np.full(self.data.num_breaths, np.nan)
        exhale_flow_peak_f = np.full(self.data.num_breaths, np.nan)

        inhale_volume_peak_t = np.full(self.data.num_breaths, np.nan)
        inhale_volume_peak_v = np.full(self.data.num_breaths, np.nan)

        exhale_volume_peak_t = np.full(self.data.num_breaths, np.nan)
        exhale_volume_peak_v = np.full(self.data.num_breaths, np.nan)

        # Setup arrays for inhale, exhale, and pause segments
        t_inhale = np.full(self.data.len_data, np.nan)
        v_inhale = np.full(self.data.len_data, np.nan)
        f_inhale = np.full(self.data.len_data, np.nan)

        t_inhale_pause = np.full(self.data.len_data, np.nan)
        v_inhale_pause = np.full(self.data.len_data, np.nan)
        f_inhale_pause = np.full(self.data.len_data, np.nan)

        t_exhale = np.full(self.data.len_data, np.nan)
        v_exhale = np.full(self.data.len_data, np.nan)
        f_exhale = np.full(self.data.len_data, np.nan)

        t_exhale_pause = np.full(self.data.len_data, np.nan)
        v_exhale_pause = np.full(self.data.len_data, np.nan)
        f_exhale_pause = np.full(self.data.len_data, np.nan)

        for idx in range(self.data.num_breaths):
            inhale_onset = self.data.inhale_onsets[idx]
            inhale_offset = self.data.inhale_offsets[idx]

            if not inhale_onset:
                # Exclude breath from being displayed
                continue

            inhale_pause_onset = self.data.inhale_pause_onsets[idx]
            inhale_pause_offset = self.data.inhale_pause_offsets[idx]

            exhale_onset = self.data.exhale_onsets[idx]
            exhale_offset = self.data.exhale_offsets[idx]

            exhale_pause_onset = self.data.exhale_pause_onsets[idx]
            exhale_pause_offset = self.data.exhale_pause_offsets[idx]

            inhale_flow_peak_t[idx] = self.data.time[self.data.inhale_flow_peaks[idx]]
            inhale_flow_peak_f[idx] = self.data.flow[self.data.inhale_flow_peaks[idx]]

            exhale_flow_peak_t[idx] = self.data.time[self.data.exhale_flow_peaks[idx]]
            exhale_flow_peak_f[idx] = self.data.flow[self.data.exhale_flow_peaks[idx]]

            inhale_volume_peak_t[idx] = self.data.time[self.data.inhale_volume_peaks[idx]]
            inhale_volume_peak_v[idx] = self.data.local_corrected_volume[self.data.inhale_volume_peaks[idx]]

            exhale_volume_peak_t[idx] = self.data.time[self.data.exhale_volume_peaks[idx]]
            exhale_volume_peak_v[idx] = self.data.local_corrected_volume[self.data.exhale_volume_peaks[idx]]

            t_inhale[inhale_onset:inhale_offset] = self.data.time[inhale_onset:inhale_offset]
            v_inhale[inhale_onset:inhale_offset] = self.data.local_corrected_volume[inhale_onset:inhale_offset]
            f_inhale[inhale_onset:inhale_offset] = self.data.flow[inhale_onset:inhale_offset]

            t_inhale_pause[inhale_pause_onset:inhale_pause_offset] = self.data.time[
                                                                     inhale_pause_onset:inhale_pause_offset]
            v_inhale_pause[inhale_pause_onset:inhale_pause_offset] = self.data.local_corrected_volume[
                                                                     inhale_pause_onset:inhale_pause_offset]
            f_inhale_pause[inhale_pause_onset:inhale_pause_offset] = self.data.flow[
                                                                     inhale_pause_onset:inhale_pause_offset]

            t_exhale[exhale_onset:exhale_offset] = self.data.time[exhale_onset:exhale_offset]
            v_exhale[exhale_onset:exhale_offset] = self.data.local_corrected_volume[exhale_onset:exhale_offset]
            f_exhale[exhale_onset:exhale_offset] = self.data.flow[exhale_onset:exhale_offset]

            t_exhale_pause[exhale_pause_onset:exhale_pause_offset] = self.data.time[
                                                                     exhale_pause_onset:exhale_pause_offset]
            v_exhale_pause[exhale_pause_onset:exhale_pause_offset] = self.data.local_corrected_volume[
                                                                     exhale_pause_onset:exhale_pause_offset]
            f_exhale_pause[exhale_pause_onset:exhale_pause_offset] = self.data.flow[
                                                                     exhale_pause_onset:exhale_pause_offset]

        # Plot volume and flow extrema
        self.plot_markers(self.ax4, inhale_volume_peak_t, inhale_volume_peak_v, False, color='k', marker=7)
        self.plot_markers(self.ax4, exhale_volume_peak_t, exhale_volume_peak_v, False, color='k', marker=6)

        self.plot_markers(self.ax5, inhale_flow_peak_t, inhale_flow_peak_f, False, color='k', marker=7)
        self.plot_markers(self.ax5, exhale_flow_peak_t, exhale_flow_peak_f, False, color='k', marker=6)

        # Plot inhale, exhale, and pause segments. Downsampling rate (dr): plot 1 in every "dr" datapoints

        self.plot_signal(self.ax4, t_inhale, v_inhale, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='b', linestyle='-')
        self.plot_signal(self.ax4, t_inhale_pause, v_inhale_pause, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='r', linestyle='-')
        self.plot_signal(self.ax4, t_exhale, v_exhale, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='g', linestyle='-')
        self.plot_signal(self.ax4, t_exhale_pause, v_exhale_pause, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='r', linestyle='-')

        self.plot_signal(self.ax5, t_inhale, f_inhale, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='b', linestyle='-')
        self.plot_signal(self.ax5, t_inhale_pause, f_inhale_pause, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='r', linestyle='-')
        self.plot_signal(self.ax5, t_exhale, f_exhale, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='g', linestyle='-')
        self.plot_signal(self.ax5, t_exhale_pause, f_exhale_pause, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=False, color='r', linestyle='-')

    def edit_calibration_point(self, event):
        self.plot_selection.edit_calibration_point(event)

    @staticmethod
    def find_closest_idx(arr, target, to_right=False, to_left=False):
        n = len(arr)
        left = 0
        right = n - 1
        mid = 0

        # edge cases:
        if target <= arr[left]:
            return left
        if target >= arr[right]:
            return right

        # find which half of array contains target, repeat
        while left < right:
            mid = (left + right) // 2
            if target < arr[mid]:
                right = mid
            elif target > arr[mid]:
                left = mid + 1
            else:
                return mid

        # if target doesn't match any array values exactly, find closest
        if target < arr[mid]:
            if to_right:
                return mid
            elif to_left:
                return mid - 1
            avg = (arr[mid - 1] + arr[mid]) / 2.0
            if target <= avg:
                return mid - 1
            else:
                return mid

        elif target > arr[mid]:
            if to_right:
                return mid + 1
            elif to_left:
                return mid
            avg = (arr[mid] + arr[mid + 1]) / 2.0
            if target <= avg:
                return mid
            else:
                return mid + 1

    def find_breath_landmarks(self):
        # ask user to define region of interest with start and stop times
        if not self.set_region_of_interest():
            return
        print('Starting Search for Breath Landmarks...')
        self.roi_start_idx, self.roi_stop_idx = self.find_start_stop_idxs(self.data.time, self.roi_start, self.roi_stop)

        # Find extrema (peaks and troughs)
        self.data.find_extrema(self.data.smoothed_flow, self.roi_start_idx, self.roi_stop_idx,
                               self.peak_detection_window_sizes, self.peak_detection_shift_sizes)

        # Find remaining breath landmarks
        self.data.find_landmarks(self.data.flow, self.data.local_corrected_volume, self.pause_search_threshold,
                                 self.number_histogram_bins, self.extra_zero_bin_threshold,
                                 self.extra_pause_bin_threshold, self.mode_bin_ratio, self.extra_pause_bin_ratio)

        # Calculate breath stats
        self.data.calculate_breath_stats()

        # Find saturated, movement, and outlier breaths
        self.data.find_saturated_breaths(self.sat_min_voltage, self.sat_max_voltage, self.remove_saturated.get())
        self.data.find_movement_breaths(self.movement_threshold, self.remove_movement.get())
        self.data.find_outlier_breaths(self.outlier_test_stats, self.outlier_window_size, self.max_num_stdevs,
                                       self.remove_outlier.get())
        self.data.determine_breaths_to_exclude(self.exclusion_adjacency, self.must_be_between_exclusions,
                                               self.remove_saturated.get(), self.remove_movement.get(),
                                               self.remove_outlier.get())

        # Calculate averages over data windows
        self.data.calculate_window_avgs(self.avg_window_size, self.avg_window_offset,
                                        self.roi_start_idx, self.roi_stop_idx)

        # Display landmarks
        self.display_breaths()

        # Notify user that search for breath landmarks has finished
        print('Search for Breath Landmarks Completed Successfully')

        # Update menu options
        self.update_menu_bar_state(self.label_edit_menu, self.label_peak_detection_window_sizes, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_peak_detection_shift_sizes, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_pause_search_threshold, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_number_histogram_bins, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_extra_zero_bin_threshold, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_extra_pause_bin_threshold, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_mode_bin_ratio, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_extra_pause_bin_ratio, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_movement_threshold, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_exclusion_adjacency, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_must_be_between_exlusions, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_avg_window_size, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_avg_window_offset, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_graphical_downsampling_rate, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_outlier_settings, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_remove_saturated, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_remove_movement, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_remove_outlier, tk.DISABLED)
        self.update_menu_bar_state(self.label_analysis_menu, self.label_find_breath_landmarks, tk.DISABLED)
        self.update_menu_bar_state(self.label_file_menu, self.label_save, tk.NORMAL)

    def find_start_stop_idxs(self, arr, start, stop):
        start_idx = self.find_closest_idx(arr, start, to_right=True)
        stop_idx = self.find_closest_idx(arr, stop, to_left=True)

        return start_idx, stop_idx

    def go_to_time(self):
        val = simpledialog.askfloat(self.label_go_to_time, 'Time in Seconds')
        if val not in self.data.time:
            return
        self.spos.val = val
        self.spos.set_val(val)

    @staticmethod
    def hide_axes(fig, ax):
        fig.delaxes(ax)

    def open(self):
        # Ask user for file to open
        new_filename = filedialog.askopenfilename()
        if new_filename == '':
            return
        # self.reinitialize()
        self.filename = new_filename
        self.short_filename = self.filename[self.filename.rfind('/') + 1:self.filename.rfind('.')]
        self.update_window_title()

        # Read data and convert voltage signal to pressure signal
        self.data = respiration.Respiration()
        self.data.load_data_file(self.filename, max_data_points=0)
        self.data.calculate_pressure_signal(self.supply_voltage, self.min_pressure, self.max_pressure)

        self.start_time = self.data.time[0]
        self.stop_time = self.data.time[-1]
        self.ax2_ymin = respiration.voltage_to_pressure(self.min_voltage, self.supply_voltage,
                                                        self.min_pressure, self.max_pressure)
        self.ax2_ymax = respiration.voltage_to_pressure(self.max_voltage, self.supply_voltage,
                                                        self.min_pressure, self.max_pressure)
        self.ax2.set_ylim(self.ax2_ymin, self.ax2_ymax)

        # Automatically calculate graphical downsampling rate. For every million datapoints above first million,
        # downsample by 1 additional datapoint, rounding up. Don't use new value if previously set value is higher.
        new_downsampling_rate = int(np.ceil(self.data.len_data / 1000000))
        if new_downsampling_rate > self.graphical_downsampling_rate:
            self.graphical_downsampling_rate = new_downsampling_rate

        # Plot data
        self.plot_signal(self.ax1, self.data.time, self.data.voltage, self.graphical_downsampling_rate,
                         replace_existing=True, show_comments=True, color='k', linestyle='-')
        self.plot_signal(self.ax2, self.data.time, self.data.pressure, self.graphical_downsampling_rate,
                         replace_existing=True, show_comments=True, color='k', linestyle='-')

        # Update menu bar options
        self.update_menu_bar_state(self.label_file_menu, self.label_open, tk.DISABLED)
        self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.NORMAL)
        self.update_menu_bar_state(self.label_analysis_menu, self.label_start_auto_calibration, tk.NORMAL)
        self.update_menu_bar_state(self.label_view_menu, self.label_go_to_time, tk.NORMAL)
        self.update_menu_bar_state(self.label_view_menu, self.label_time_scale, tk.NORMAL)
        self.update_menu_bar_state(self.label_view_menu, self.label_voltage_scale, tk.NORMAL)
        self.update_menu_bar_state(self.label_view_menu, self.label_pressure_scale, tk.NORMAL)

    def plot_signal(self, ax, x, y, downsampling_rate, replace_existing=True, show_comments=True, **kwargs):
        if replace_existing:
            for line in ax.get_lines():
                line.remove()
        new_line, = ax.plot(x[::downsampling_rate], y[::downsampling_rate], **kwargs)
        if show_comments:
            for data_idx in self.data.comment_idxs:
                ax.annotate(self.data.comments[data_idx], (x[data_idx], y[data_idx]), xycoords='data', xytext=(-30, 50),
                            textcoords='offset pixels', horizontalalignment='right', verticalalignment='top',
                            arrowprops=dict(arrowstyle='-|>', facecolor='k',
                                            connectionstyle='angle,angleA=0,angleB=90'))

        self.spos.valmin = self.start_time
        self.spos.valmax = self.stop_time - self.time_window
        self.spos.ax.set_xlim(self.spos.valmin, self.spos.valmax)
        self.canvas.draw()
        return new_line

    def plot_markers(self, ax, x, y, replace_existing=False, **kwargs):
        if replace_existing:
            for line in ax.get_lines():
                line.remove()
        path_collection = ax.scatter(x, y, **kwargs)

        self.spos.valmin = self.start_time
        self.spos.valmax = self.stop_time - self.time_window
        self.spos.ax.set_xlim(self.spos.valmin, self.spos.valmax)
        self.canvas.draw()
        return path_collection

    def save(self):
        filename = filedialog.asksaveasfilename(confirmoverwrite=True,
                                                initialdir=self.filename[0:self.filename.rfind('/') + 1],
                                                initialfile=self.short_filename + '_analysis',
                                                defaultextension='.xlsx',
                                                filetypes=[('Excel Workbook', '*.xlsx')])
        if filename:
            print("Saving Analysis...")
            self.data.write_data_to_xlsx(filename, self)
            short_filename = filename[filename.rfind('/') + 1:filename.rfind('.')]
            messagebox.showinfo('File Saved', '%s.xlsx saved successfully.' % short_filename)
            print('%s.xlsx Saved Successfully' % short_filename)

    def set_autocal_settings(self):
        AutocalibrationSettingsDialog(self, title='Autocalibration Settings')

    def set_avg_window_offset(self):
        value = simpledialog.askfloat('Set Averaging Window Offset', 'Offset [s]:',
                                      initialvalue=self.avg_window_offset, minvalue=0.0)
        if value is not None:
            self.avg_window_offset = value

    def set_avg_window_size(self):
        value = simpledialog.askfloat('Set Averaging Window Size', 'Window Size [s]:',
                                      initialvalue=self.avg_window_size, minvalue=0.0)
        if value is not None:
            self.avg_window_size = value

    def set_downsampling_rate(self):
        value = simpledialog.askinteger('Set Graphical Downsampling Rate', 'Downsampling Rate:',
                                        initialvalue=self.graphical_downsampling_rate, minvalue=1)

        if value is not None:
            self.graphical_downsampling_rate = value

    def set_exclusion_adjacency(self):
        value = simpledialog.askinteger('Set Breath Exclusion Adjacency', 'Number of Breaths:',
                                        initialvalue=self.exclusion_adjacency, minvalue=0)
        if value is not None:
            self.exclusion_adjacency = value

    def set_extra_pause_bin_ratio(self):
        value = simpledialog.askfloat('Set Pause Search Extra Pause Bin Ratio', 'Ratio:',
                                      initialvalue=self.extra_pause_bin_ratio, minvalue=0.0, maxvalue=1.0)
        if value is not None:
            self.extra_pause_bin_ratio = value

    def set_extra_pause_bin_threshold(self):
        value = simpledialog.askfloat('Set Pause Search Extra Pause Bin Threshold', 'Threshold [%]:',
                                      initialvalue=self.extra_pause_bin_threshold * 100.0, minvalue=0.0, maxvalue=100.0)
        if value is not None:
            self.extra_pause_bin_threshold = value / 100.0

    def set_extra_zero_bin_threshold(self):
        value = simpledialog.askfloat('Set Pause Search Extra Zero Bin Threshold', 'Threshold [%]:',
                                      initialvalue=self.extra_zero_bin_threshold * 100.0, minvalue=0.0, maxvalue=100.0)
        if value is not None:
            self.extra_zero_bin_threshold = value / 100.0

    def set_filter_cutoff(self):
        value = simpledialog.askfloat('Set Lowpass Filter Cutoff Frequency', 'Lowpass Cutoff Frequency [Hz]:',
                                      initialvalue=self.volume_calibration_cutoff, minvalue=0.0)
        if value is not None:
            self.volume_calibration_cutoff = value

    def set_filter_order(self):
        value = simpledialog.askinteger('Set Lowpass Filter Order', 'Lowpass Filter Order:',
                                        initialvalue=self.volume_calibration_order, minvalue=1)
        if value is not None:
            self.volume_calibration_order = value

    def set_flow_oversmoothing_sigma(self):
        value = simpledialog.askfloat('Set Flow Signal Gaussian Smoothing Filter Sigma',
                                      'Gaussian Smoothing Filter Sigma:',
                                      initialvalue=self.flow_oversmoothing_sigma, minvalue=0.0)
        if value is not None:
            self.flow_oversmoothing_sigma = value

    def set_local_drift_corr_window_size(self):
        value = simpledialog.askfloat('Set Local Drift Correction Window Size', 'Window Size [s]:',
                                      initialvalue=self.local_drift_correction_window_size, minvalue=0.0)
        if value is not None:
            self.local_drift_correction_window_size = value

    def set_mode_bin_ratio(self):
        value = simpledialog.askfloat('Set Pause Search Mode Bin Ratio', 'Ratio:',
                                      initialvalue=self.mode_bin_ratio, minvalue=1.0)
        if value is not None:
            self.mode_bin_ratio = value

    def set_movement_threshold(self):
        value = simpledialog.askfloat('Set Movement Threshold', 'Movement Threshold [mL]:',
                                      initialvalue=self.movement_threshold, minvalue=0.0)
        if value is not None:
            self.movement_threshold = value

    def set_must_be_between_exclusions(self):
        if self.must_be_between_exclusions:
            initial_value = 1
        else:
            initial_value = 0
        value = simpledialog.askinteger('Set Double-Sided Exclusion Adjacency Requirement', '0 - False, 1 - True:',
                                        initialvalue=initial_value, minvalue=0, maxvalue=1)
        if value is not None:
            if value == 1:
                self.must_be_between_exclusions = True
            elif value == 0:
                self.must_be_between_exclusions = False

    def set_number_histogram_bins(self):
        value = simpledialog.askinteger('Set Pause Search Number of Histogram Bins', 'Number of Bins:',
                                        initialvalue=self.number_histogram_bins, minvalue=0)
        if value is not None:
            self.number_histogram_bins = value

    def set_outlier_settings(self):
        OutlierSettingsDialog(self, title='Outlier Settings')

    def set_pause_search_threshold(self):
        value = simpledialog.askfloat('Set Pause Search Volume Threshold', 'Threshold [%]:',
                                      initialvalue=self.pause_search_threshold, minvalue=0.0, maxvalue=100.0)
        if value is not None:
            self.pause_search_threshold = value

    def set_peak_detection_shift_sizes(self):
        shift_sizes = ','.join([str(format(x * 100.0, '.1f')) for x in self.peak_detection_shift_sizes])
        values = simpledialog.askstring('Set Peak Detection Shift Sizes', 'Shift Sizes [%]:',
                                        initialvalue=shift_sizes)
        if values is not None:
            split_values = values.split(',')
            new_window_sizes = [float(x.strip()) / 100.0 for x in split_values]
            self.peak_detection_shift_sizes = new_window_sizes

    def set_peak_detection_window_sizes(self):
        window_sizes = ','.join([str(x) for x in self.peak_detection_window_sizes])
        values = simpledialog.askstring('Set Peak Detection Window Sizes', 'Window Sizes [samples]:',
                                        initialvalue=window_sizes)
        if values is not None:
            split_values = values.split(',')
            new_window_sizes = [int(x.strip()) for x in split_values]
            self.peak_detection_window_sizes = new_window_sizes

    def set_region_of_interest(self):
        value1 = simpledialog.askfloat('Set Region of Interest', 'Start Time [s]:',
                                       initialvalue=self.data.time[0], minvalue=0.0)
        if value1 is not None:
            self.roi_start = value1
        else:
            return False

        value2 = simpledialog.askfloat('Set Region of Interest', 'Stop Time [s]:',
                                       initialvalue=self.data.time[-1], minvalue=value1)

        if value2 is not None:
            self.roi_stop = value2
        else:
            return False

        self.roi_start_idx, self.roi_stop_idx = self.find_start_stop_idxs(self.data.time, self.roi_start, self.roi_stop)
        return True

    def set_sat_max_voltage(self):
        value = simpledialog.askfloat('Set Upper Saturation Voltage', 'Upper Saturation Voltage [V]:',
                                      initialvalue=self.sat_max_voltage, minvalue=self.min_voltage,
                                      maxvalue=self.max_voltage)
        if value is not None:
            self.sat_max_voltage = value

    def set_sat_min_voltage(self):
        value = simpledialog.askfloat('Set Lower Saturation Voltage', 'Lower Saturation Voltage [V]:',
                                      initialvalue=self.sat_min_voltage, minvalue=self.min_voltage,
                                      maxvalue=self.max_voltage)
        if value is not None:
            self.sat_min_voltage = value

    def set_smoothing_sigma(self):
        value = simpledialog.askfloat('Set Volume Signal Gaussian Smoothing Filter Sigma',
                                      'Gaussian Smoothing Filter Sigma:',
                                      initialvalue=self.volume_smoothing_sigma, minvalue=0.0)
        if value is not None:
            self.volume_smoothing_sigma = value

    @staticmethod
    def show_axes(fig, ax):
        fig.add_subplot(ax)

    @staticmethod
    def smooth_flow_signal(flow_signal, smoothing_sigma):
        smoothed = respiration.smoothing(flow_signal, smoothing_sigma)
        return smoothed

    def start_auto_calibration(self):
        # Verify that a pressure signal exists.
        if not self.data.pressure.size:
            return

        # Ask user if they want to use calibration data from a previous file
        self.use_prev_cal = messagebox.askyesno('Volume Calibration', 'Use calibration data from another file?')
        if self.use_prev_cal:
            self.prev_cal_filename = filedialog.askopenfilename()
            if self.prev_cal_filename != '':
                print("Using Previous Calibration Data...")
                min_pressure, max_pressure, ref_pressure = self.data.read_calibration_datafile(self.prev_cal_filename)
                self.min_pressure = min_pressure
                self.max_pressure = max_pressure
                self.ref_pressure = ref_pressure

                # Recalculate pressure signal to use min/max/ref pressures from previous calibration file
                self.data.calculate_pressure_signal(self.supply_voltage, self.min_pressure, self.max_pressure)

                # Update menu bar options
                self.update_menu_bar_state(self.label_view_menu, self.label_voltage_scale, tk.DISABLED)
                self.update_menu_bar_state(self.label_view_menu, self.label_volume_scale, tk.NORMAL)
                self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.DISABLED)
                self.update_menu_bar_state(self.label_analysis_menu, self.label_start_auto_calibration, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_sat_min_voltage, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_sat_max_voltage, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_filter_cutoff, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_filter_order, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_autocal_settings, tk.DISABLED)

                # Finish calibration automatically
                self.stop_volume_calibration()

            else:
                return

        # Otherwise, let program auto-calibrate
        else:
            print("Starting Auto-Calibration...")
            self.data.lowpass_filter_pressure(self.volume_calibration_cutoff, self.volume_calibration_order)

            self.hide_axes(self.fig, self.ax1)
            self.show_axes(self.fig, self.ax3)
            self.fig.canvas.draw_idle()

            # Use MSF algorithm to find steps in roi's around calibration injections,
            # should we use raw or lowpass filtered pressure signal?
            self.data.moving_step_fit(self.data.pressure, self.autocal_sigma, self.autocal_search_left,
                                      self.autocal_search_right, self.autocal_dup_comment_distance,
                                      self.autocal_lowest_cal_vol, self.autocal_highest_cal_vol,
                                      self.autocal_ref_cal_vol, self.autocal_w_ref_cal_vol, self.autocal_w_variance,
                                      self.autocal_num_std_step_height_threshold, self.autocal_step_width_threshold,
                                      self.autocal_lower_variance_ratio, self.autocal_downsampling_rate,
                                      self.autocal_num_std_exclude, self.autocal_make_plots)

            # Plot calibration points and fit line
            path1 = self.plot_markers(self.ax3, self.data.excluded_calibration_volumes,
                                      self.data.excluded_calibration_pressures, replace_existing=True, color='r')
            self.plot_signal(self.ax3, [self.ax3_xmin, self.ax3_xmax],
                             [self.data.excluded_calibration_m * self.ax3_xmin,
                              self.data.excluded_calibration_m * self.ax3_xmax],
                             downsampling_rate=1, replace_existing=False,
                             show_comments=False, color='r', linestyle='--')
            path2 = self.plot_markers(self.ax3, self.data.calibration_volumes, self.data.calibration_pressures,
                                      replace_existing=False, color='k')
            line2 = self.plot_signal(self.ax3, [self.ax3_xmin, self.ax3_xmax],
                                     [self.data.calibration_m * self.ax3_xmin,
                                      self.data.calibration_m * self.ax3_xmax],
                                     downsampling_rate=1, replace_existing=False,
                                     show_comments=False, color='k', linestyle='-')
            path1.set_label('%s Excluded Points' % len(self.data.excluded_calibration_volumes))
            path2.set_label('%s Included Points' % len(self.data.calibration_volumes))
            line2.set_label('y = %.6f x\nr^2 = %.6f' % (self.data.calibration_m, self.data.calibration_r_squared))
            self.ax3.legend(loc='upper left')

            self.canvas.draw()
            self.check_calibration_axes()

            # Update menu bar options
            self.update_menu_bar_state(self.label_view_menu, self.label_voltage_scale, tk.DISABLED)
            self.update_menu_bar_state(self.label_view_menu, self.label_inj_volume_scale, tk.NORMAL)
            self.update_menu_bar_state(self.label_view_menu, self.label_cal_pressure_scale, tk.NORMAL)
            self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.DISABLED)
            self.update_menu_bar_state(self.label_analysis_menu, self.label_start_auto_calibration, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_sat_min_voltage, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_sat_max_voltage, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_filter_cutoff, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_filter_order, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_autocal_settings, tk.DISABLED)

            # Finish calibration automatically
            self.stop_volume_calibration()

    def start_volume_calibration(self):
        # Verify that a pressure signal exists.
        if not self.data.pressure.size:
            return

        # Ask user if they want to use calibration data from a previous file
        self.use_prev_cal = messagebox.askyesno('Volume Calibration', 'Use calibration data from another file?')

        if self.use_prev_cal:
            self.prev_cal_filename = filedialog.askopenfilename()
            if self.prev_cal_filename != '':
                print("Using Previous Calibration Data...")
                min_pressure, max_pressure, ref_pressure = self.data.read_calibration_datafile(self.prev_cal_filename)
                self.min_pressure = min_pressure
                self.max_pressure = max_pressure
                self.ref_pressure = ref_pressure

                # Recalculate pressure signal to use min/max/ref pressures from previous calibration file
                self.data.calculate_pressure_signal(self.supply_voltage, self.min_pressure, self.max_pressure)

                # Update menu bar options
                self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.DISABLED)

                # Update menu bar options
                self.update_menu_bar_state(self.label_view_menu, self.label_voltage_scale, tk.DISABLED)
                self.update_menu_bar_state(self.label_view_menu, self.label_volume_scale, tk.NORMAL)
                self.update_menu_bar_state(self.label_analysis_menu, self.label_start_auto_calibration, tk.DISABLED)
                self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.DISABLED)
                self.update_menu_bar_state(self.label_analysis_menu, self.label_stop_volume_calibration, tk.NORMAL)
                self.update_menu_bar_state(self.label_edit_menu, self.label_sat_min_voltage, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_sat_max_voltage, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_filter_cutoff, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_filter_order, tk.DISABLED)
                self.update_menu_bar_state(self.label_edit_menu, self.label_autocal_settings, tk.DISABLED)

                # Finish calibration automatically
                self.stop_volume_calibration()

            else:
                return

        else:
            print("Starting Manual Calibration...")
            # Filter and plot pressure signal with a lowpass Butterworth filter
            self.data.lowpass_filter_pressure(self.volume_calibration_cutoff, self.volume_calibration_order)
            line = self.plot_signal(self.ax2, self.data.time, self.data.lowpass_pressure,
                                    self.graphical_downsampling_rate, replace_existing=False,
                                    show_comments=False, color='r', linestyle='-')
            self.hide_axes(self.fig, self.ax1)
            self.show_axes(self.fig, self.ax3)
            self.fig.canvas.draw_idle()

            # Allow user to select calibration peaks
            self.plot_selection = PlotSelection(self, self.data, line)

            # Update menu bar options
            self.update_menu_bar_state(self.label_view_menu, self.label_voltage_scale, tk.DISABLED)
            self.update_menu_bar_state(self.label_view_menu, self.label_inj_volume_scale, tk.NORMAL)
            self.update_menu_bar_state(self.label_view_menu, self.label_cal_pressure_scale, tk.NORMAL)
            self.update_menu_bar_state(self.label_analysis_menu, self.label_start_auto_calibration, tk.DISABLED)
            self.update_menu_bar_state(self.label_analysis_menu, self.label_start_volume_calibration, tk.DISABLED)
            self.update_menu_bar_state(self.label_analysis_menu, self.label_stop_volume_calibration, tk.NORMAL)
            self.update_menu_bar_state(self.label_edit_menu, self.label_sat_min_voltage, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_sat_max_voltage, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_filter_cutoff, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_filter_order, tk.DISABLED)
            self.update_menu_bar_state(self.label_edit_menu, self.label_autocal_settings, tk.DISABLED)

    def stop_volume_calibration(self):
        if self.use_prev_cal:
            done = True
        else:
            done = messagebox.askyesno('Stop Volume Calibration', 'Finish and Save Volume Calibration?')
        if not done:
            return

        # Update saturation values
        self.sat_min_pressure = respiration.voltage_to_pressure(self.sat_min_voltage, self.supply_voltage,
                                                                self.min_pressure, self.max_pressure)
        self.sat_max_pressure = respiration.voltage_to_pressure(self.sat_max_voltage, self.supply_voltage,
                                                                self.min_pressure, self.max_pressure)
        self.sat_min_volume = respiration.pressure_to_volume(self.sat_min_pressure, self.ref_pressure,
                                                             self.data.calibration_m, self.data.calibration_b)
        self.sat_max_volume = respiration.pressure_to_volume(self.sat_max_pressure, self.ref_pressure,
                                                             self.data.calibration_m, self.data.calibration_b)

        # Calculate volume signal using calibration data
        self.data.calculate_volume_signal(self.ref_pressure)

        # Smooth volume signal using Gaussian filter
        self.data.smooth_volume_signal(self.volume_smoothing_sigma)

        # Todo: Give user option to toggle global and local drift correction
        # Remove global drift using linear regression on whole data set
        self.data.remove_global_drift(self.data.time, self.data.smoothed_volume)

        # Remove local drift using mean filter with 2 second window
        self.data.remove_local_drift(self.data.global_corrected_volume, self.local_drift_correction_window_size)

        # Plot volume signal
        if self.use_prev_cal:
            self.hide_axes(self.fig, self.ax1)
        else:
            self.hide_axes(self.fig, self.ax3)
        self.show_axes(self.fig, self.ax4)
        self.plot_signal(self.ax4, self.data.time, self.data.local_corrected_volume, self.graphical_downsampling_rate,
                         replace_existing=False, show_comments=True, color='k', linestyle='-')

        # Reset min/max volume range based on min/max pressure values
        self.axis_min_pressure = respiration.voltage_to_pressure(self.min_voltage, self.supply_voltage,
                                                                 self.min_pressure, self.max_pressure)
        self.axis_max_pressure = respiration.voltage_to_pressure(self.max_voltage, self.supply_voltage,
                                                                 self.min_pressure, self.max_pressure)
        self.min_volume = respiration.pressure_to_volume(self.axis_min_pressure, self.ref_pressure,
                                                         self.data.calibration_m, self.data.calibration_b)
        self.max_volume = respiration.pressure_to_volume(self.axis_max_pressure, self.ref_pressure,
                                                         self.data.calibration_m, self.data.calibration_b)

        self.ax4_ymin = self.min_volume
        self.ax4_ymax = self.max_volume
        self.ax4.set_ylim(self.ax4_ymin, self.ax4_ymax)
        self.fig.canvas.draw_idle()

        # Calculate flow signal
        self.calculate_volumetric_flow()

        # Oversmooth flow signal to more easily detect peaks in flow signal
        smoothed_flow = self.smooth_flow_signal(self.data.flow, self.flow_oversmoothing_sigma)
        self.data.smoothed_flow = smoothed_flow
        if self.display_smoothed_flow.get():
            self.plot_signal(self.ax5, self.data.time, smoothed_flow, self.graphical_downsampling_rate,
                             replace_existing=False, show_comments=False, color='c', linestyle='-')
            self.fig.canvas.draw_idle()

        # Notify user that calibration completed successfully
        print("Calibration Completed Successfully")

        # Update menu bar options
        self.update_menu_bar_state(self.label_edit_menu, self.label_volume_smoothing_sigma, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_local_drift_corr_window_size, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_flow_oversmoothing_sigma, tk.DISABLED)
        self.update_menu_bar_state(self.label_edit_menu, self.label_display_smoothed_flow, tk.DISABLED)
        self.update_menu_bar_state(self.label_analysis_menu, self.label_stop_volume_calibration, tk.DISABLED)
        self.update_menu_bar_state(self.label_view_menu, self.label_inj_volume_scale, tk.DISABLED)
        self.update_menu_bar_state(self.label_view_menu, self.label_cal_pressure_scale, tk.DISABLED)
        self.update_menu_bar_state(self.label_view_menu, self.label_volume_scale, tk.NORMAL)

    def update_menu_bar_state(self, toolbar_menu, menu_to_update, new_state):
        index = self.submenus[toolbar_menu].index(menu_to_update)
        self.submenus[toolbar_menu].entryconfigure(index, state=new_state)

    def update_window_title(self):
        if not self.short_filename:
            self.winfo_toplevel().title(self.title)
        else:
            self.winfo_toplevel().title('%s - %s' % (self.short_filename, self.title))


class AutocalibrationSettingsDialog(simpledialog.Dialog):
    def __init__(self, parent, title=None):
        simpledialog.Dialog.__init__(self, parent, title)
        # the following initialized variables will be further initialized in the body method
        self.var_smoothing_sigma = None
        self.var_search_left = None
        self.var_search_right = None
        self.var_dup_comment_distance = None
        self.var_lowest_cal_vol = None
        self.var_highest_cal_vol = None
        self.var_ref_cal_vol = None
        self.var_w_ref_cal_vol = None
        self.var_w_variance = None
        self.var_num_std_step_height_threshold = None
        self.var_step_width_threshold = None
        self.var_lower_variance_ratio = None
        self.var_downsampling_rate = None
        self.var_num_std_exclude = None
        self.cb_make_plots_state = None
        self.cb_make_plots = None

        self.entry_smoothing_sigma = None
        self.entry_search_left = None
        self.entry_search_right = None
        self.entry_dup_comment_distance = None
        self.entry_lowest_cal_vol = None
        self.entry_highest_cal_vol = None
        self.entry_ref_cal_vol = None
        self.entry_w_ref_cal_vol = None
        self.entry_w_variance = None
        self.entry_num_std_step_height_threshold = None
        self.entry_step_width_threshold = None
        self.entry_lower_variance_ratio = None
        self.entry_downsampling_rate = None
        self.entry_num_std_exclude = None

    def body(self, parent):
        label_smoothing_sigma = tk.Label(parent, text='Gaussian Filter Smoothing Sigma [s]')
        label_search_left = tk.Label(parent, text='Seconds to Search Left of Comment [s]')
        label_search_right = tk.Label(parent, text='Seconds to Search Right of Comment [s]')
        label_dup_comment_distance = tk.Label(parent, text='Duplicate Comment Distance [s]')
        label_lowest_cal_vol = tk.Label(parent, text='Lowest Calibration Volume [mL]')
        label_highest_cal_vol = tk.Label(parent, text='Highest Calibration Volume [mL]')
        label_ref_cal_vol = tk.Label(parent, text='Reference Calibration Volume [mL]')
        label_w_ref_cal_vol = tk.Label(parent, text='Step Window Width for Reference Calibration Volume [s]')
        label_w_variance = tk.Label(parent, text='Variance Window Width [s]')
        label_num_std_step_height_threshold = tk.Label(parent,
                                                       text='Number of Standard Deviations for Step Height Threshold')
        label_step_width_threshold = tk.Label(parent, text='Step Width Threshold [s]')
        label_lower_variance_ratio = tk.Label(parent, text='Lower Variance Ratio')
        label_downsampling_rate = tk.Label(parent, text='Downsampling Rate')
        label_num_std_exclude = tk.Label(parent, text='Number of Standard Deviations for Calibration Point Exclusion')

        self.var_smoothing_sigma = tk.StringVar()
        self.var_search_left = tk.StringVar()
        self.var_search_right = tk.StringVar()
        self.var_dup_comment_distance = tk.StringVar()
        self.var_lowest_cal_vol = tk.StringVar()
        self.var_highest_cal_vol = tk.StringVar()
        self.var_ref_cal_vol = tk.StringVar()
        self.var_w_ref_cal_vol = tk.StringVar()
        self.var_w_variance = tk.StringVar()
        self.var_num_std_step_height_threshold = tk.StringVar()
        self.var_step_width_threshold = tk.StringVar()
        self.var_lower_variance_ratio = tk.StringVar()
        self.var_downsampling_rate = tk.StringVar()
        self.var_num_std_exclude = tk.StringVar()
        self.cb_make_plots_state = tk.IntVar()

        self.entry_smoothing_sigma = tk.Entry(parent, textvariable=self.var_smoothing_sigma)
        self.entry_search_left = tk.Entry(parent, textvariable=self.var_search_left)
        self.entry_search_right = tk.Entry(parent, textvariable=self.var_search_right)
        self.entry_dup_comment_distance = tk.Entry(parent, textvariable=self.var_dup_comment_distance)
        self.entry_lowest_cal_vol = tk.Entry(parent, textvariable=self.var_lowest_cal_vol)
        self.entry_highest_cal_vol = tk.Entry(parent, textvariable=self.var_highest_cal_vol)
        self.entry_ref_cal_vol = tk.Entry(parent, textvariable=self.var_ref_cal_vol)
        self.entry_w_ref_cal_vol = tk.Entry(parent, textvariable=self.var_w_ref_cal_vol)
        self.entry_w_variance = tk.Entry(parent, textvariable=self.var_w_variance)
        self.entry_num_std_step_height_threshold = tk.Entry(parent, textvariable=self.var_num_std_step_height_threshold)
        self.entry_step_width_threshold = tk.Entry(parent, textvariable=self.var_step_width_threshold)
        self.entry_lower_variance_ratio = tk.Entry(parent, textvariable=self.var_lower_variance_ratio)
        self.entry_downsampling_rate = tk.Entry(parent, textvariable=self.var_downsampling_rate)
        self.entry_num_std_exclude = tk.Entry(parent, textvariable=self.var_num_std_exclude)

        self.var_smoothing_sigma.set(str(self.parent.autocal_sigma))
        self.var_search_left.set(str(self.parent.autocal_search_left))
        self.var_search_right.set(str(self.parent.autocal_search_right))
        self.var_dup_comment_distance.set(str(self.parent.autocal_dup_comment_distance))
        self.var_lowest_cal_vol.set(str(self.parent.autocal_lowest_cal_vol))
        self.var_highest_cal_vol.set(str(self.parent.autocal_highest_cal_vol))
        self.var_ref_cal_vol.set(str(self.parent.autocal_ref_cal_vol))
        self.var_w_ref_cal_vol.set(str(self.parent.autocal_w_ref_cal_vol))
        self.var_w_variance.set(str(self.parent.autocal_w_variance))
        self.var_num_std_step_height_threshold.set(str(self.parent.autocal_num_std_step_height_threshold))
        self.var_step_width_threshold.set(str(self.parent.autocal_step_width_threshold))
        self.var_lower_variance_ratio.set(str(self.parent.autocal_lower_variance_ratio))
        self.var_downsampling_rate.set(str(self.parent.autocal_downsampling_rate))
        self.var_num_std_exclude.set(str(self.parent.autocal_num_std_exclude))
        self.cb_make_plots = tk.Checkbutton(parent, text='Show All Autocalibration Plots and Data',
                                            variable=self.cb_make_plots_state)
        if self.parent.autocal_make_plots:
            self.cb_make_plots.select()
        else:
            self.cb_make_plots.deselect()

        label_smoothing_sigma.grid(sticky=tk.W, padx=5)
        self.entry_smoothing_sigma.grid(sticky=tk.W + tk.E, padx=5)
        label_search_left.grid(sticky=tk.W, padx=5)
        self.entry_search_left.grid(sticky=tk.W + tk.E, padx=5)
        label_search_right.grid(sticky=tk.W, padx=5)
        self.entry_search_right.grid(sticky=tk.W + tk.E, padx=5)
        label_dup_comment_distance.grid(sticky=tk.W, padx=5)
        self.entry_dup_comment_distance.grid(sticky=tk.W + tk.E, padx=5)
        label_lowest_cal_vol.grid(sticky=tk.W, padx=5)
        self.entry_lowest_cal_vol.grid(sticky=tk.W + tk.E, padx=5)
        label_highest_cal_vol.grid(sticky=tk.W, padx=5)
        self.entry_highest_cal_vol.grid(sticky=tk.W + tk.E, padx=5)
        label_ref_cal_vol.grid(sticky=tk.W, padx=5)
        self.entry_ref_cal_vol.grid(sticky=tk.W + tk.E, padx=5)
        label_w_ref_cal_vol.grid(sticky=tk.W, padx=5)
        self.entry_w_ref_cal_vol.grid(sticky=tk.W + tk.E, padx=5)
        label_w_variance.grid(sticky=tk.W, padx=5)
        self.entry_w_variance.grid(sticky=tk.W + tk.E, padx=5)
        label_num_std_step_height_threshold.grid(sticky=tk.W, padx=5)
        self.entry_num_std_step_height_threshold.grid(sticky=tk.W + tk.E, padx=5)
        label_step_width_threshold.grid(sticky=tk.W, padx=5)
        self.entry_step_width_threshold.grid(sticky=tk.W + tk.E, padx=5)
        label_lower_variance_ratio.grid(sticky=tk.W, padx=5)
        self.entry_lower_variance_ratio.grid(sticky=tk.W + tk.E, padx=5)
        label_downsampling_rate.grid(sticky=tk.W, padx=5)
        self.entry_downsampling_rate.grid(sticky=tk.W + tk.E, padx=5)
        label_num_std_exclude.grid(sticky=tk.W, padx=5)
        self.entry_num_std_exclude.grid(sticky=tk.W + tk.E, padx=5)
        self.cb_make_plots.grid(sticky=tk.W, padx=15)

        return self.entry_smoothing_sigma

    def apply(self):
        self.parent.autocal_sigma = float(self.entry_smoothing_sigma.get())
        self.parent.autocal_search_left = float(self.entry_search_left.get())
        self.parent.autocal_search_right = float(self.entry_search_right.get())
        self.parent.autocal_dup_comment_distance = float(self.entry_dup_comment_distance.get())
        self.parent.autocal_lowest_cal_vol = float(self.entry_lowest_cal_vol.get())
        self.parent.autocal_highest_cal_vol = float(self.entry_highest_cal_vol.get())
        self.parent.autocal_ref_cal_vol = float(self.entry_ref_cal_vol.get())
        self.parent.autocal_w_ref_cal_vol = float(self.entry_w_ref_cal_vol.get())
        self.parent.autocal_w_variance = float(self.entry_w_variance.get())
        self.parent.autocal_num_std_step_height_threshold = float(self.entry_num_std_step_height_threshold.get())
        self.parent.autocal_step_width_threshold = float(self.entry_step_width_threshold.get())
        self.parent.autocal_lower_variance_ratio = float(self.entry_lower_variance_ratio.get())
        self.parent.autocal_downsampling_rate = int(self.entry_downsampling_rate.get())
        self.parent.autocal_num_std_exclude = float(self.entry_num_std_exclude.get())
        self.parent.autocal_make_plots = bool(self.cb_make_plots_state.get())


class OutlierSettingsDialog(simpledialog.Dialog):
    def __init__(self, parent, title=None):
        simpledialog.Dialog.__init__(self, parent, title)
        # the following initialized variables will be further initialized in the body method
        self.cb_piv_state = None
        self.cb_pev_state = None
        self.cb_pif_state = None
        self.cb_pef_state = None
        self.cb_ti_state = None
        self.cb_te_state = None
        self.cb_tb_state = None
        self.entry_outlier_window_size = None
        self.entry_outlier_threshold = None
        self.cb_piv = None
        self.cb_pev = None
        self.cb_pif = None
        self.cb_pef = None
        self.cb_ti = None
        self.cb_te = None
        self.cb_tb = None
        self.entry1 = None
        self.entry2 = None

    def body(self, parent):
        self.cb_piv_state = tk.IntVar()
        self.cb_pev_state = tk.IntVar()
        self.cb_pif_state = tk.IntVar()
        self.cb_pef_state = tk.IntVar()
        self.cb_ti_state = tk.IntVar()
        self.cb_te_state = tk.IntVar()
        self.cb_tb_state = tk.IntVar()
        self.entry_outlier_window_size = tk.StringVar()
        self.entry_outlier_threshold = tk.StringVar()
        label1 = tk.Label(parent, text='Parameters to Examine for Outliers:')
        label2_text = 'Outlier Determination Window Size [s]:'
        label2 = tk.Label(parent, text=label2_text)
        label3_text = 'Maximum Number of Standard Deviations from Mean:'
        label3 = tk.Label(parent, text=label3_text)
        self.cb_piv = tk.Checkbutton(parent, text='Peak Inspiratory Volume', variable=self.cb_piv_state)
        self.cb_pev = tk.Checkbutton(parent, text='Peak Expiratory Volume', variable=self.cb_pev_state)
        self.cb_pif = tk.Checkbutton(parent, text='Peak Inspiratory Flow', variable=self.cb_pif_state)
        self.cb_pef = tk.Checkbutton(parent, text='Peak Expiratory Flow', variable=self.cb_pef_state)
        self.cb_ti = tk.Checkbutton(parent, text='Duration of Inspiration', variable=self.cb_ti_state)
        self.cb_te = tk.Checkbutton(parent, text='Duration of Expiration', variable=self.cb_te_state)
        self.cb_tb = tk.Checkbutton(parent, text='Duration of Breath', variable=self.cb_tb_state)
        self.entry1 = tk.Entry(parent, textvariable=self.entry_outlier_window_size)
        outlier_window_size = str(self.parent.outlier_window_size)
        self.entry_outlier_window_size.set(outlier_window_size)
        self.entry2 = tk.Entry(parent, textvariable=self.entry_outlier_threshold)
        outlier_threshold = str(self.parent.max_num_stdevs)
        self.entry_outlier_threshold.set(outlier_threshold)
        if 'piv' in self.parent.outlier_test_stats:
            self.cb_piv.select()
        else:
            self.cb_piv.deselect()
        if 'pev' in self.parent.outlier_test_stats:
            self.cb_pev.select()
        else:
            self.cb_pev.deselect()
        if 'pif' in self.parent.outlier_test_stats:
            self.cb_pif.select()
        else:
            self.cb_pif.deselect()
        if 'pef' in self.parent.outlier_test_stats:
            self.cb_pef.select()
        else:
            self.cb_pef.deselect()
        if 'ti' in self.parent.outlier_test_stats:
            self.cb_ti.select()
        else:
            self.cb_ti.deselect()
        if 'te' in self.parent.outlier_test_stats:
            self.cb_te.select()
        else:
            self.cb_te.deselect()
        if 'tb' in self.parent.outlier_test_stats:
            self.cb_tb.select()
        else:
            self.cb_tb.deselect()

        label1.grid(sticky=tk.W, padx=5)
        self.cb_piv.grid(sticky=tk.W, padx=15)
        self.cb_pev.grid(sticky=tk.W, padx=15)
        self.cb_pif.grid(sticky=tk.W, padx=15)
        self.cb_pef.grid(sticky=tk.W, padx=15)
        self.cb_ti.grid(sticky=tk.W, padx=15)
        self.cb_te.grid(sticky=tk.W, padx=15)
        self.cb_tb.grid(sticky=tk.W, padx=15)
        label2.grid(sticky=tk.W, padx=5)
        self.entry1.grid(sticky=tk.W + tk.E, padx=5)
        label3.grid(sticky=tk.W, padx=5)
        self.entry2.grid(sticky=tk.W + tk.E, padx=5)
        return self.cb_piv

    def apply(self):
        outlier_test_stats = []
        if bool(self.cb_piv_state.get()):
            outlier_test_stats.append('piv')
        if bool(self.cb_pev_state.get()):
            outlier_test_stats.append('pev')
        if bool(self.cb_pif_state.get()):
            outlier_test_stats.append('pif')
        if bool(self.cb_pef_state.get()):
            outlier_test_stats.append('pef')
        if bool(self.cb_ti_state.get()):
            outlier_test_stats.append('ti')
        if bool(self.cb_te_state.get()):
            outlier_test_stats.append('te')
        if bool(self.cb_tb_state.get()):
            outlier_test_stats.append('tb')

        self.parent.outlier_test_stats = outlier_test_stats
        self.parent.outlier_window_size = float(self.entry_outlier_window_size.get())
        self.parent.max_num_stdevs = float(self.entry_outlier_threshold.get())


class PlotSelection:
    def __init__(self, parent, respiration_data, line):
        self.parent = parent
        self.data = respiration_data
        self.line = line
        self.sampling_rate = respiration_data.sampling_rate
        self.ax = self.line.axes
        self.fig = self.ax.get_figure()
        self.line_selections = []
        self.xdata = self.line.get_xdata()
        self.ydata = self.line.get_ydata()
        self.idxs = []
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.ann = None
        self.ann_coord = tuple()

    def __call__(self, event):
        if event.xdata is None or self.xdata is None:
            return False
        if event.inaxes != self.line.axes:
            return False
        idx = self.find_nearest_idx(event.xdata)
        self.add_idx(idx)

    def add_idx(self, idx):
        # user selects 4 points in total, two for the zero portion, and two for the injection portion of the calibration
        # each potion is averaged and the difference between the averages of these portions is reported
        idxs = list(self.idxs)
        if len(idxs) >= 4:
            idxs = []
            self.clear_points()
        elif idx in idxs:
            return
        idxs.append(idx)
        self.idxs = idxs
        self.draw_selection_point(idx)
        if len(self.idxs) == 4:
            idxs.sort()
            self.find_calibration_point()

    def clear_points(self):
        for line in self.line_selections:
            line.remove()
        self.line_selections = []
        self.fig.canvas.draw_idle()

    def draw_selection_point(self, idx):
        (x, y) = self.get_point_from_idx(idx)
        line, = self.ax.plot(x, y, 'rx')

        self.line_selections.append(line)
        self.fig.canvas.draw_idle()

    def edit_calibration_point(self, event):
        if event.artist == self.parent.ax3_line1:
            line = event.artist
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            ind = event.ind

            # Highlight selected calibration point
            selection_line, = self.parent.ax3.plot(np.take(xdata, ind)[0], np.take(ydata, ind)[0],
                                                   'ro', fillstyle='none', markersize=10)
            self.fig.canvas.draw_idle()

            # Ask user if they want to remove selected calibration point
            m = messagebox.askyesno('Edit Calibration Points', 'Remove selected calibration point?')
            selection_line.remove()
            if m:
                removal_coord = tuple([self.data.calibration_volumes[int(ind[0])],
                                       self.data.calibration_pressures[int(ind[0])]])
                self.data.remove_calibration_point(int(ind[0]))

                self.parent.ax3_line1.set_data(self.data.calibration_volumes, self.data.calibration_pressures)
                if self.ann and self.ann_coord == removal_coord:
                    self.ann.remove()
                    self.ann = None
                    self.ann_coord = tuple()

                # There should be two or more points, with two or more unique x values to perform the polyfit
                if len(set(self.data.calibration_volumes)) >= 2:
                    self.polyfit_calibration_data()
                else:
                    self.data.calibration_m = None
                    self.data.calibration_b = None
                    self.parent.ax3_line2.set_data([], [])
                    legend = self.parent.ax3.get_legend()
                    if legend:
                        legend.remove()

            self.fig.canvas.draw_idle()

    def find_calibration_point(self):
        if len(self.idxs) != 4:
            return
        idx1 = self.idxs[0]
        idx2 = self.idxs[1]
        idx3 = self.idxs[2]
        idx4 = self.idxs[3]

        # compute averages
        p1_avg = np.average(self.ydata[idx1:idx2])
        p2_avg = np.average(self.ydata[idx3:idx4])

        p1_raw_avg = np.average(self.data.pressure[idx1:idx2])
        p2_raw_avg = np.average(self.data.pressure[idx3:idx4])

        # Ask user for injection volume
        volume = simpledialog.askfloat('Calibration Injection Point', 'Injection Volume in mLs:')
        if not volume:
            # Clear selections and extrema point
            self.clear_points()
            return
        else:
            # Verify the sign of the volume change is correct.
            # For positive injections: dV_chamber < 0, dV_lungs = -dV_chamber, dV_lungs > 0, dP > 0.
            # For negative injections: dV_chamber > 0, dV_lungs = -dV_chamber, dV_lungs < 0, dP < 0.
            if p2_raw_avg >= p1_raw_avg:
                volume = abs(volume)
            else:
                volume = (-1.0) * abs(volume)

        calibration_time = self.xdata[idx3]
        calibration_pressure = p2_avg - p1_avg
        calibration_volume = volume

        # Troubleshooting
        # p_max_raw = max(self.data.pressure[idx1:idx4])

        # Add calibration point to Respiration data object
        self.data.add_calibration_point(calibration_time, calibration_pressure, calibration_volume)

        # Add calibration point to calibration axis of Analysis Window
        self.parent.ax3_line1.set_data(self.data.calibration_volumes, self.data.calibration_pressures)
        coord_str = '(%.2f, %.0f)' % (calibration_volume, calibration_pressure)
        if self.ann:
            self.ann.remove()
        self.ann_coord = tuple([calibration_volume, calibration_pressure])
        self.ann = self.parent.ax3.annotate(coord_str, (calibration_volume, calibration_pressure),
                                            xycoords='data', xytext=(0, 20),
                                            textcoords='offset pixels', horizontalalignment='center',
                                            verticalalignment='top')
        self.parent.check_calibration_axes()

        # There should be two or more points, with two or more unique x values to perform the polyfit
        if len(set(self.data.calibration_volumes)) >= 2:
            self.polyfit_calibration_data()

        self.fig.canvas.draw_idle()

    def find_nearest_idx(self, value):
        idx = np.searchsorted(self.xdata, value, side="left")
        if idx > 0 and (idx == len(self.xdata) or abs(value - self.xdata[idx - 1]) < abs(value - self.xdata[idx])):
            return idx - 1
        else:
            return idx

    def get_point_from_idx(self, idx):
        time = self.xdata[idx]
        pressure = self.ydata[idx]
        return tuple([time, pressure])

    def polyfit_calibration_data_old(self):
        z = np.polyfit(self.data.calibration_volumes, self.data.calibration_pressures, 1)
        m, b = float(z[0]), float(z[1])
        p = np.poly1d(z)
        # Calculate r_squared value
        y_mean = np.average(self.data.calibration_pressures)
        unexplained_error = 0
        total_error = 0
        for idx in range(len(self.data.calibration_volumes)):
            y = self.data.calibration_pressures[idx]
            y_pred = p(self.data.calibration_volumes[idx])
            unexplained_error += (y - y_pred) ** 2
            total_error += (y - y_mean) ** 2
        r_squared = 1 - (unexplained_error / total_error)
        self.parent.ax3_line2.set_data([self.parent.ax3_xmin, self.parent.ax3_xmax],
                                       [p(self.parent.ax3_xmin), p(self.parent.ax3_xmax)])
        self.parent.ax3_line2.set_label('y = %.6f x + %.6f\nr^2 = %.6f' % (m, b, r_squared))
        self.parent.ax3.legend(loc='upper left')
        self.data.calibration_m = m
        self.data.calibration_b = b
        self.data.calibration_r_squared = r_squared
        return

    def polyfit_calibration_data(self):
        # Fit calibration points with line through origin (y = mx).
        # m = sum(x_i * y_i) / sum(x_i * x_i)
        times = self.data.calibration_times
        volumes = self.data.calibration_volumes
        pressures = self.data.calibration_pressures
        numerator_sum = 0
        denominator_sum = 0
        for time, volume, pressure in zip(times, volumes, pressures):
            numerator_i = volume * pressure
            denominator_i = volume ** 2
            numerator_sum += numerator_i
            denominator_sum += denominator_i
        if denominator_sum == 0:
            slope = 0
        else:
            slope = numerator_sum / denominator_sum

        # Calculate r_squared value
        mean_pressure = np.mean(pressures)
        unexplained_error = 0
        total_error = 0
        for volume, pressure in zip(volumes, pressures):
            expected_pressure = slope * volume
            unexplained_error += (pressure - expected_pressure) ** 2
            total_error += (pressure - mean_pressure) ** 2
        r_squared = 1 - (unexplained_error / total_error)

        # Plot regression line and equation with r_squared value
        self.parent.ax3_line2.set_data([self.parent.ax3_xmin, self.parent.ax3_xmax],
                                       [slope * self.parent.ax3_xmin, slope * self.parent.ax3_xmax])
        self.parent.ax3_line2.set_label('y = %.6f x\nr^2 = %.6f' % (slope, r_squared))
        self.parent.ax3.legend(loc='upper left')

        # Save regression parameters
        self.data.calibration_m = slope
        self.data.calibration_b = 0
        self.data.calibration_r_squared = r_squared

        return


def main():
    root = tk.Tk()
    analysis_window = AnalysisWindow(root)
    root.protocol('WM_DELETE_WINDOW', analysis_window.close_window)
    root.mainloop()


if __name__ == '__main__':
    main()
