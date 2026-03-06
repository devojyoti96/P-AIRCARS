import numpy as np
from casatools import table, msmetadata
from joblib import Parallel, delayed
import time
from .imaging import calc_sun_dia

   
def uvsub_flagger(uvwave, data, threshold=5.0, min_bin_samples=5, num_bins=50, binning_type="log"):
    """
    Flag data based on binning in uv-wavelengths space, specifically designed for solar data.
    Optimized with Numba.

    Args:
        uvwave        : UV distances in wavelengths (1D numpy array)
        data          : Visibility amplitudes (1D numpy array for this optimized version)
        threshold     : Flagging threshold multiplier for MAD (default: 5.0)
        min_bin_samples: Minimum number of samples required in a bin (default: 5)
        num_bins      : Number of logarithmic bins to use (default: 50)

    Returns:
        Boolean array with same shape as data indicating flagged points
    """
    if uvwave.shape[0] == 0 or data.shape[0] == 0:
        return np.zeros(data.shape, dtype=np.bool_)

    if uvwave.shape[0] != data.shape[0]:
        return np.zeros(data.shape, dtype=np.bool_)

    nrows = data.shape[0]

    flags = np.zeros(nrows, dtype=np.bool_)

    # Handle case with too few valid UV points or data points
    valid_mask = ~np.isnan(data) & ~np.isnan(uvwave)
    n_valid = np.sum(valid_mask)
    if n_valid < min_bin_samples * 2:
        print("Returning")
        return flags 

    valid_uvwave = uvwave[valid_mask]
    valid_data = data[valid_mask]

    # Create non-uniform bins
    positive_uv = valid_uvwave[valid_uvwave > 0]
    if len(positive_uv) == 0:
        return flags 

    min_uv = np.min(positive_uv)
    max_uv = np.max(valid_uvwave) 

    # Check for invalid range before logspace
    if min_uv <= 0 or max_uv <= min_uv:
        return flags 

    # Adaptive binning: more bins where data density is higher 
    current_num_bins = num_bins
    if max_uv / min_uv > 1000:
        current_num_bins = min(100, num_bins * 2)

    # Ensure num_bins is at least 2 for logspace
    if current_num_bins < 2:
        current_num_bins = 2

    # bins = np.logspace(np.log10(min_uv), np.log10(max_uv), current_num_bins)
    if binning_type == "log":
        log_bins = np.linspace(np.log10(min_uv), np.log10(max_uv), current_num_bins)
        bins = 10**log_bins
    else:
        bins = np.linspace(min_uv, max_uv, current_num_bins)

    # Bin the data and compute statistics
    # Note: np.digitize behavior with Numba might need verification, especially edge cases.
    # Indices correspond to bins[i-1] <= x < bins[i]
    # We want bins based on bin edges, so adjust indices.
    bin_indices = np.searchsorted(bins, valid_uvwave, side='right')

    # For each bin, compute median and MAD
    # Use numba.prange for potential parallelization within the function if needed,
    # but primary parallelism is outside over timestamps.
    for bin_idx in range(current_num_bins): # Iterate through bin indices 0 to num_bins-1
        # Find data points belonging to this bin
        # bin_mask corresponds to points where bin_indices == bin_idx
        bin_mask = (bin_indices == bin_idx)
        bin_count = np.sum(bin_mask)

        if bin_count < min_bin_samples:  # Skip bins with too few points
            continue

        bin_data = valid_data[bin_mask]

        # Robust statistics
        # Use trimmed statistics if enough points 
        if bin_count > 20:
            # Sort data and trim extreme values
            sorted_data = np.sort(bin_data)
            trim_size = max(1, int(bin_count * 0.05))  # Trim 5% from each end
            # Ensure trim_size doesn't exceed half the array size
            if 2 * trim_size >= bin_count:
                 trim_size = max(0, (bin_count - 1) // 2) # Adjust trim size to avoid empty array

            if trim_size > 0 :
                 trimmed_data = sorted_data[trim_size:-trim_size]
                 if trimmed_data.size == 0: # Check if trimming resulted in empty array
                     bin_median = np.median(bin_data) # Fallback to full median
                 else:
                     bin_median = np.median(trimmed_data)
            else:
                 bin_median = np.median(bin_data)

        else:
            bin_median = np.median(bin_data)

        # Use robust MAD calculation
        abs_diff = np.abs(bin_data - bin_median)
        bin_mad = np.median(abs_diff)

        # Avoid division by zero or very small MAD
        mad_threshold = 1e-9 # Adjusted threshold for comparison
        if bin_mad < mad_threshold:
            data_range = np.max(bin_data) - np.min(bin_data)
            if data_range < 1e-6: # Adjusted threshold
                # Very uniform data, skip flagging
                continue
            else:
                # Use a fraction of the data range as MAD
                bin_mad = data_range * 0.1
                if bin_mad < mad_threshold: # Ensure MAD is not too small even after adjustment
                    bin_mad = mad_threshold

        # Adaptive threshold based on bin size and data spread
        adaptive_threshold = threshold
        if bin_count < 10:
            # Be more conservative with small bins
            adaptive_threshold *= 1.5

        # Flag outliers in this bin
        # Use the calculated adaptive_threshold and bin_mad
        # Ensure bin_mad is positive before division
        if bin_mad > 0:
             deviation = abs_diff / bin_mad
             outliers_mask = deviation > adaptive_threshold
        else:
             # Handle case where bin_mad is effectively zero after checks
             # Only flag if deviation from median is large in absolute terms?
             # Or skip flagging in this case. Let's skip for now.
             outliers_mask = np.zeros(bin_data.shape, dtype=np.bool_)


        # Only flag if not too many points would be flagged (avoid over-flagging)
        outlier_fraction = np.sum(outliers_mask) / bin_count if bin_count > 0 else 0.0

        final_outliers_mask = outliers_mask # Start with the initial outlier mask

        if outlier_fraction > 0.4:  # If more than 40% would be flagged
            # This might be a real feature, not RFI - be more conservative
            if bin_mad > 0: # Check bin_mad again
                 extreme_outliers_mask = deviation > (adaptive_threshold * 2)
                 if np.sum(extreme_outliers_mask) / bin_count < 0.2:  # If less than 20% are extreme outliers
                     final_outliers_mask = extreme_outliers_mask # Use the more conservative mask
                 else:
                     # Skip flagging this bin - likely a real feature
                     final_outliers_mask = np.zeros(bin_data.shape, dtype=np.bool_) # Clear flags for this bin
            else:
                 # Skip flagging if MAD is zero and outlier fraction is high
                 final_outliers_mask = np.zeros(bin_data.shape, dtype=np.bool_)


        # Map back to original indices within the valid data subset
        # Find indices within valid_data that correspond to this bin
        bin_indices_in_valid_data = np.where(bin_mask)[0]
        # Apply the final outlier mask to these indices
        outlier_indices_in_bin = bin_indices_in_valid_data[final_outliers_mask]

        # Map these indices back to the original data array size (nrows)
        # Find the original indices corresponding to the valid data
        original_indices_of_valid_data = np.where(valid_mask)[0]
        # Get the original indices of the outliers
        original_outlier_indices = original_indices_of_valid_data[outlier_indices_in_bin]

        # Set flags in the main flags array
        flags[original_outlier_indices] = True

    return flags


# Helper function for parallel processing of a single timestamp
def _process_timestamp(time_indices, uv_distances_full, data_masked_full, nchan, npol, threshold=3.0, num_bins=25, binning_type='log'):
    """Processes a single timestamp's data to generate flags."""
    timestamp_flags = np.zeros((len(time_indices), nchan, npol), dtype=bool)
    timestamp_uv = uv_distances_full[time_indices]

    for chan in range(nchan):
        for pol in range(npol):
            # Get data for this timestamp, channel, and polarization
            # Extract data only for the relevant indices
            timestamp_data = np.abs(data_masked_full[time_indices, chan, pol])

            # Skip if too many NaNs or too few points
            valid_data_mask = ~np.isnan(timestamp_data)
            n_valid = np.sum(valid_data_mask)
            if n_valid < 10: # Use a reasonable minimum number of points
                continue

            # Apply UV-based flagging only on this timestamp's valid data
            # Pass only the valid subset to uvsub_flagger
            current_flags = uvsub_flagger(
                timestamp_uv[valid_data_mask],
                timestamp_data[valid_data_mask],
                threshold=threshold,
                min_bin_samples=5, # Keep min_bin_samples consistent
                num_bins=num_bins,
                binning_type=binning_type 
            )

            # Update flags for this timestamp, channel, and polarization
            if current_flags is not None and len(current_flags) > 0:
                 # Map flags back to the original size of the timestamp slice
                 valid_indices_in_timestamp = np.where(valid_data_mask)[0]
                 timestamp_flags[valid_indices_in_timestamp[current_flags], chan, pol] = True

    return time_indices, timestamp_flags


def flagger(msname, datacolumn, threshold=3.0, num_processes=4, num_bins=30, binning_type="log", visualize=False): # Added num_bins
    """
    Applies the hybrid solar model flagging approach to the measurement set.
    Optimized for speed using Numba and parallel processing.

    Args:
        msname       : Path to the measurement set.
        datacolumn   : Name of the data column (e.g. 'DATA' or 'CORRECTED_DATA').
        threshold    : Multiplier for the MAD-based flagging threshold.
        num_processes: Number of processes for parallel processing.
        num_bins     : Number of UV bins for uvsub_flagger (default: 30)
    """
    print(f"Flagging : {msname}")
    import time
    start_time = time.time()
    ms = table()
    msmd = msmetadata()
    msmd.open(msname)
    freq_Hz = msmd.meanfreq(0)
    msmd.close()
    ms.open(msname, nomodify=False) #readonly=False, ack=False)
    try:
        time_col = ms.getcol("TIME")
        unique_times = np.unique(time_col)
        #print(f"Number of unique timestamps: {len(unique_times)}")

        # --- Determine Data Shape ---
        # Read the data column first to get its actual shape
        #print(f"Reading data column '{datacolumn}' to determine shape...")
        if datacolumn == "RESIDUAL":
            data = ms.getcol("CORRECTED_DATA") - ms.getcol("MODEL_DATA")
        else:
            data = ms.getcol(datacolumn)
        data = data.T
        #print("Data reading complete.")

        data_actual_shape = data.shape
        #print(f"Actual data shape: {data_actual_shape}")

        if len(data_actual_shape) == 3:
            n_rows, nchan, npol = data_actual_shape
            #print(f"Determined shape: rows={n_rows}, nchan={nchan}, npol={npol}")
        else:
            # Handle unexpected dimensions
            raise ValueError(
                f"Unexpected data dimensions in column '{datacolumn}'. "
                f"Expected 3 (rows, chans, pols), got {len(data_actual_shape)} with shape {data_actual_shape}."
            )

        # --- Get or Create FLAG Column ---
        if "FLAG" in ms.colnames():
            #print("Reading existing FLAG column...")
            flags = ms.getcol("FLAG").T
            # Check if flag shape matches data shape
            if flags.shape != data_actual_shape:
                raise ValueError(
                    f"FLAG column shape {flags.shape} does not match data column shape {data_actual_shape}."
                )
        else:
            #print("Creating new FLAG array (in memory)...")
            flags = np.zeros(data_actual_shape, dtype=bool) # Use determined shape

        # --- Get UVW Data ---
        uvw = ms.getcol("UVW").T
        #print (uvw.shape)
        if uvw.shape[0] != n_rows:
             raise ValueError(f"UVW table rows ({uvw.shape[0]}) do not match data rows ({n_rows})")

        wavelength = 299792458 / freq_Hz
        uvw_wavelength = uvw / wavelength


        # --- Initial Flag Count ---
        # Number of unflagged data points initially
        n_flagged = np.sum(flags)
        n_total = flags.size
        n_unflagged = n_total - n_flagged
        #print(f"Initial flag count: {n_flagged}/{n_total} ({n_flagged/n_total:.2%})")

        data[flags]=np.nan
        
        # Calculate UV distances in wavelengths
        uv_distances = np.sqrt(uvw_wavelength[:, 0]**2 + uvw_wavelength[:, 1]**2)

        # --- Parallel Processing Setup ---
        #print(f"Processing timestamps in parallel using {num_processes} workers...")

        # Pre-calculate indices for each timestamp to avoid repeated searching
        # Original slow method:
        # time_indices_map = {ts: np.where(time_col == ts)[0] for ts in unique_times}

        # Faster method using argsort:
        #print("Building timestamp-to-index map (using argsort)...")
        sort_start_time = time.time()
        # Get indices that would sort the time column
        sorted_indices = np.argsort(time_col, kind='stable') # Stable sort preserves original order for ties
        # Get the sorted times and the locations where the time changes
        sorted_times = time_col[sorted_indices]
        unique_times_sorted, split_points = np.unique(sorted_times, return_index=True)
        # Split the sorted_indices array at these change points
        # This gives lists of original indices, grouped by timestamp
        indices_split_by_time = np.split(sorted_indices, split_points[1:])
        # Create the map from unique time to its corresponding array of indices
        time_indices_map = dict(zip(unique_times_sorted, indices_split_by_time))
        sort_end_time = time.time()
        #print(f"Timestamp map built in {sort_end_time - sort_start_time:.2f} seconds.")

        # Verify unique_times consistency (optional check)
        if not np.array_equal(np.sort(unique_times), unique_times_sorted):
            print("Warning: Unique times from initial scan and sorting differ. Using sorted unique times.")
            unique_times = unique_times_sorted # Use the times derived from sorting


        # Prepare arguments for each parallel job
        job_args = [
            (
                time_indices_map[timestamp], 
                uv_distances,                
                data,                        
                nchan,
                npol,
                threshold,
                num_bins,
                binning_type                    
            )
            for timestamp in unique_times if len(time_indices_map[timestamp]) >= 10 # Filter small timestamps here
        ]

        #print('Done with job_args')

        #print(f"Number of timestamps to process in parallel: {len(job_args)}")
        if not job_args:
            print("Warning: No timestamps met the minimum data requirement for parallel processing.")
            # Handle case with no jobs gracefully, maybe skip parallel section
            results = [] 
        elif len(job_args) < num_processes:
            print(f"Warning: Number of jobs ({len(job_args)}) is less than requested processes ({num_processes}). Effective parallelism will be limited.")

        #print('Done with job_args check')

        parallel_start_time = time.time()
        if job_args: # Only run if there are jobs
            results = Parallel(n_jobs=num_processes, backend='loky')( 
                delayed(_process_timestamp)(*args) for args in job_args
            )
        parallel_end_time = time.time()
        #print(f"Parallel processing part took: {parallel_end_time - parallel_start_time:.2f} seconds")


        # --- Combine Results ---
        combine_start_time = time.time()
        #print("Combining results from parallel workers...")
        # Create a copy of flags to update, or update in place if acceptable
        new_flags = flags.copy()
        for time_indices, timestamp_flags in results:
            # Update the main flags array using the indices and the returned flags
            # Ensure shapes match before assignment
            if new_flags[time_indices, :, :].shape == timestamp_flags.shape:
                 new_flags[time_indices, :, :] = np.logical_or(new_flags[time_indices, :, :], timestamp_flags)
            else:
                 print(f"Warning: Shape mismatch when combining flags for indices {time_indices}. Skipping update.")
                 print(f"Expected shape: {new_flags[time_indices, :, :].shape}, Got shape: {timestamp_flags.shape}")

        combine_end_time = time.time()
        #print(f"Combining results took: {combine_end_time - combine_start_time:.2f} seconds")

        # Number of additional flagged data points
        n_final_flagged = np.sum(new_flags)
        n_additional_flagged = n_final_flagged - n_flagged
        #print(f"Number of additional flagged data points: {n_additional_flagged}")

        # Fraction of initially unflagged data points that are now flagged
        flag_fraction = n_additional_flagged / n_unflagged if n_unflagged > 0 else 0
        #print(f"Fraction of initially unflagged data flagged: {flag_fraction:.4f}")
        #print(f"Total flagged fraction: {n_final_flagged / n_total:.4f}")


        end_time = time.time()
        #print(f"Flagging logic completed in: {end_time - start_time:.2f} seconds")

        # --- Visualization (Optional) ---
        if visualize:
            print("Generating visualization plots...")
            from matplotlib import pyplot as plt
            plt.switch_backend('TkAgg')

        if visualize:
            # Plot the flagged data for a specific channel/pol
            debug_chan = min(6, nchan - 1) 
            debug_pol = min(0, npol - 1)   

            #fig = plt.figure(figsize=(16, 8))
            #ax1 = fig.add_subplot(2, 2, 1)
            #ax2 = fig.add_subplot(2, 2, 2)
            #ax3 = fig.add_subplot(2, 2, 3)
            #ax4 = fig.add_subplot(2, 2, 4)

            original_data_abs = np.abs(data[:, debug_chan, debug_pol])
            original_data_phase = np.angle(data[:, debug_chan, debug_pol]) * 180 / np.pi

            initially_unflagged_mask = ~flags[:, debug_chan, debug_pol]

            finally_unflagged_mask = ~new_flags[:, debug_chan, debug_pol]

            '''ax1.set_title(f'Initially Unflagged Amplitude (Ch {debug_chan}, Pol {debug_pol})')
            ax1.scatter(uv_distances[initially_unflagged_mask],
                        original_data_abs[initially_unflagged_mask],
                        marker='.', color='blue', alpha=0.3, s=5) # Smaller points
            ax1.set_xlabel('UV Distance (λ)')
            ax1.set_ylabel('Amplitude')
            ax1.grid(True, which='both', linestyle='--', alpha=0.5)
            ax1.set_yscale('linear') # Use log scale for amplitude often

            # Second subplot: Initially unflagged data phase
            ax2.set_title(f'Initially Unflagged Phase (Ch {debug_chan}, Pol {debug_pol})')
            ax2.scatter(uv_distances[initially_unflagged_mask],
                        original_data_phase[initially_unflagged_mask],
                        marker='.', color='blue', alpha=0.3, s=5)
            ax2.set_xlabel('UV Distance (λ)')
            ax2.set_ylabel('Phase (deg)')
            ax2.grid(True, which='both', linestyle='--', alpha=0.5)
            ax2.set_ylim(-180, 180) # Keep phase range consistent

            # Third subplot: Finally unflagged data amplitude
            ax3.set_title(f'Finally Unflagged Amplitude (Ch {debug_chan}, Pol {debug_pol})')
            ax3.scatter(uv_distances[finally_unflagged_mask],
                        original_data_abs[finally_unflagged_mask],
                        marker='.', color='green', alpha=0.3, s=5)
            ax3.set_xlabel('UV Distance (λ)')
            ax3.set_ylabel('Amplitude')
            ax3.grid(True, which='both', linestyle='--', alpha=0.5)
            ax3.set_yscale('linear') # Match scale with ax1

            # Fourth subplot: Finally unflagged data phase
            ax4.set_title(f'Finally Unflagged Phase (Ch {debug_chan}, Pol {debug_pol})')
            ax4.scatter(uv_distances[finally_unflagged_mask],
                        original_data_phase[finally_unflagged_mask],
                        marker='.', color='green', alpha=0.3, s=5)
            ax4.set_xlabel('UV Distance (λ)')
            ax4.set_ylabel('Phase (deg)')
            ax4.grid(True, which='both', linestyle='--', alpha=0.5)
            ax4.set_ylim(-180, 180) # Match scale with ax2

            plt.tight_layout()
            plot_filename = f"{msname}_flagger_vis_ch{debug_chan}_pol{debug_pol}.png"
            #plt.savefig(plot_filename) # Save plot instead of showing interactively
            #print(f"Visualization saved to {plot_filename}")
            #plt.close(fig) # Close figure to free memory
            plt.show()'''

            newly_flagged_mask = new_flags[:, debug_chan, debug_pol] & ~flags[:, debug_chan, debug_pol]
            if np.any(newly_flagged_mask):
                 fig_new = plt.figure(figsize=(10, 6))
                 plt.title(f'Newly Flagged Points (Ch {debug_chan}, Pol {debug_pol})')
                 # plot the unflagged points in blue
                 plt.scatter(uv_distances[finally_unflagged_mask], original_data_abs[finally_unflagged_mask], marker='.', color='blue', alpha=0.5, s=20, label='Unflagged')
                 # plot the newly flagged points in red
                 plt.scatter(uv_distances[newly_flagged_mask],
                             original_data_abs[newly_flagged_mask],
                             marker='.', color='red', alpha=0.5, s=20, label='Newly Flagged')
                 plt.xlabel('UV Distance (λ)')
                 plt.ylabel('Amplitude')
                 plt.xscale('linear')
                 plt.yscale('linear')
                 plt.grid(True, which='both', linestyle='--', alpha=0.5)
                 plt.legend()
                 plot_filename_new = f"{msname}_flagger_newly_flagged_ch{debug_chan}_pol{debug_pol}.png"
                 #plt.savefig(plot_filename_new)
                 #print(f"Newly flagged points plot saved to {plot_filename_new}")
                 #plt.close(fig_new)
                 plt.show()


        # Update the FLAG column in the MS
        flagput_start_time = time.time()
        #print("Updating FLAG column in MS...")
        ms.putcol("FLAG", new_flags.T)
        flagput_end_time = time.time()
        #print(f"FLAG column updated in {flagput_end_time - flagput_start_time:.2f} seconds")

    finally:
        ms.close()
        #print("MS table closed.")
    print("Solar UV-based flagging complete.")

    # Remove the second interactive plot show() call
    # plt.figure(figsize=(16, 8))
    # # Donot plot the 0 uvw points
    # plt.scatter(uv_distances[uv_distances > 0], np.abs(data[uv_distances > 0,6,0]), marker='.', color='red', alpha=0.5)
    # plt.xlabel('UV Distance (m)')
    # plt.ylabel('Amplitude')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.show()
    return 0

