import plotly.graph_objects as go
from spec import *

def plot_spectra_comparison(input_list, lwl=4000, rwl=6000, degrade=None, stack_spectra=True, vertical_offset=0.3):
    """
    Plots a list of stellar spectra together using Plotly for interactive comparison.

    Parameters
    ----------
    input_list : str or list
        Filename of the list (e.g., 'BSGs_CNO_IDs.txt') or a Python list of star IDs/names.
    lwl : float, optional
        Left wavelength limit in Angstroms. Default is 4000 Å.
    rwl : float, optional
        Right wavelength limit in Angstroms. Default is 6000 Å.
    stack_spectra : bool, optional
        If True, stacks spectra vertically. If False, overlays them on the same baseline.
    vertical_offset : float, optional
        Vertical spacing between spectra when stacked. Default is 0.3.
    """
    # 1. Resolve the input into a list of stars
    if isinstance(input_list, str):
        if '.txt' in input_list or '.lst' in input_list:
            stars = findlist(input_list)
        else:
            stars = [input_list]
    elif isinstance(input_list, list):
        stars = input_list
    else:
        print("Error: input_list must be a filename string or a list of star IDs.")
        return

    # 2. Create the Plotly figure
    fig = go.Figure()
    loaded_count = 0

    for i, star in enumerate(stars):
        try:
            # Load spectrum with the best SNR
            spectrum = spec(star, snr='best')

            # Crop to specified wavelength range
            spectrum.waveflux(lwl=lwl, rwl=rwl)

            if spectrum is None or spectrum.wave is None or spectrum.flux is None:
                print(f"Skipping {star}: No spectrum data found.")
                continue

            wave = spectrum.wave
            flux = spectrum.flux

            # Calculate plotted flux based on stacking configuration
            plotted_flux = flux + (loaded_count * vertical_offset) if stack_spectra else flux

            # Custom hover template
            hover_tpl = (
                f"<b>{star}</b><br>"
                "Wavelength: %{x:.2f} Å<br>"
                "Plotted Y: %{y:.3f}<br>"
            )
            if stack_spectra:
                hover_tpl += "True Flux: %{customdata:.3f}<extra></extra>"
            else:
                hover_tpl += "<extra></extra>"

            # Add spectrum trace
            fig.add_trace(go.Scatter(
                x=wave,
                y=plotted_flux,
                mode='lines',
                name=star,
                line=dict(width=1.2),
                customdata=flux,  # Stores original unshifted flux for tooltips
                hovertemplate=hover_tpl
            ))
            print(f"Successfully loaded {star}")
            loaded_count += 1

        except Exception as e:
            print(f"Error loading {star}: {e}")

    if loaded_count == 0:
        print("No spectra were successfully loaded. Plot will not be shown.")
        return

    # 3. Customize layout for scientific comparison
    fig.update_layout(
        title={
            'text': "Quick-View of Stellar Spectra",
            'y': 0.95,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title="Wavelength (Å)",
        yaxis_title="Normalized Flux" + (" (+ Offset)" if stack_spectra else ""),
        template="plotly_white",
        hovermode="x unified" if not stack_spectra else "closest",
        legend=dict(
            title="Stars",
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    # 4. Show the interactive plot
    fig.show()


# --- Example Usage ---
if __name__ == "__main__":
    # Example 1: Load from text file, cropped between 4000 and 4600 Angstroms
    plot_spectra_comparison(
        input_list='Benchmark_BSGs.txt',
        lwl=4000,
        rwl=4600,
        stack_spectra=True,
        vertical_offset=0.25
    )

    # Example 2: Pass a list directly and overlay them on top of each other
    # my_stars = ['HD37022', 'HD34078']
    # plot_spectra_comparison(my_stars, lwl=4500, rwl=4900, stack_spectra=False)