import os
import threading
from pathlib import Path

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError as exc:  # pragma: no cover - handled at runtime for missing Tk
    tk = None  # type: ignore[assignment]
    filedialog = None  # type: ignore[assignment]
    messagebox = None  # type: ignore[assignment]

from .pipeline import run_pipeline


def _pick_alleles(entry: tk.Entry) -> None:
    path = filedialog.askopenfilename(
        title="Select allele FASTA file",
        filetypes=[("FASTA files", "*.fa *.fna *.fasta"), ("All files", "*.*")],
    )
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


def _pick_alleles_dir(entry: tk.Entry) -> None:
    path = filedialog.askdirectory(title="Select allele FASTA directory")
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


def _pick_genome_dir(entry: tk.Entry) -> None:
    path = filedialog.askdirectory(title="Select genome FASTA directory")
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


def _pick_output(entry: tk.Entry) -> None:
    path = filedialog.asksaveasfilename(
        title="Select output CSV",
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


def main() -> int:
    if tk is None:
        raise SystemExit(
            "Tkinter is not available in this Python installation. "
            "Please install a Python build with Tk support."
        )

    root = tk.Tk()
    root.title("mbcgmlst")
    root.resizable(False, False)

    button_width = 9

    alleles_label = tk.Label(root, text="Allele FASTA (file or directory):")
    alleles_entry = tk.Entry(root, width=60)
    alleles_buttons = tk.Frame(root)
    alleles_file_btn = tk.Button(
        alleles_buttons,
        text="File...",
        width=button_width,
        command=lambda: _pick_alleles(alleles_entry),
    )
    alleles_dir_btn = tk.Button(
        alleles_buttons,
        text="Dir...",
        width=button_width,
        command=lambda: _pick_alleles_dir(alleles_entry),
    )

    genome_label = tk.Label(root, text="Genome FASTA directory:")
    genome_entry = tk.Entry(root, width=60)
    genome_btn = tk.Button(
        root, text="Browse...", width=button_width, command=lambda: _pick_genome_dir(genome_entry)
    )

    output_label = tk.Label(root, text="Output CSV:")
    output_entry = tk.Entry(root, width=60)
    output_btn = tk.Button(
        root, text="Save...", width=button_width, command=lambda: _pick_output(output_entry)
    )

    threads_label = tk.Label(root, text="Threads:")
    threads_var = tk.StringVar(value=str(os.cpu_count() or 1))
    threads_frame = tk.Frame(root)
    threads_entry = tk.Entry(threads_frame, width=10, textvariable=threads_var)
    threads_max_btn = tk.Button(
        threads_frame,
        text="Max",
        width=button_width,
        command=lambda: threads_var.set(str(os.cpu_count() or 1)),
    )

    progress_label = tk.Label(root, text="Progress:")
    progress_var = tk.IntVar(value=0)
    progress_bar = ttk.Progressbar(
        root, orient="horizontal", mode="determinate", variable=progress_var, maximum=1
    )

    status_var = tk.StringVar(value="Idle")
    status_label = tk.Label(root, textvariable=status_var, anchor="w")

    run_button = tk.Button(root, text="Run", width=12)

    root.grid_columnconfigure(1, weight=1)

    alleles_label.grid(row=0, column=0, sticky="w", padx=8, pady=(8, 6))
    alleles_entry.grid(row=0, column=1, sticky="we", padx=8, pady=(8, 6))
    alleles_buttons.grid(row=0, column=2, sticky="w", padx=4, pady=(8, 6))
    alleles_file_btn.grid(row=0, column=0, padx=(0, 4))
    alleles_dir_btn.grid(row=0, column=1)

    genome_label.grid(row=1, column=0, sticky="w", padx=8, pady=6)
    genome_entry.grid(row=1, column=1, sticky="we", padx=8, pady=6)
    genome_btn.grid(row=1, column=2, sticky="w", padx=4, pady=6)

    output_label.grid(row=2, column=0, sticky="w", padx=8, pady=6)
    output_entry.grid(row=2, column=1, sticky="we", padx=8, pady=6)
    output_btn.grid(row=2, column=2, sticky="w", padx=4, pady=6)

    threads_label.grid(row=3, column=0, sticky="w", padx=8, pady=6)
    threads_frame.grid(row=3, column=1, sticky="w", padx=8, pady=6)
    threads_entry.grid(row=0, column=0)
    threads_max_btn.grid(row=0, column=1, padx=(6, 0))

    progress_label.grid(row=4, column=0, sticky="w", padx=8, pady=6)
    progress_bar.grid(row=4, column=1, sticky="we", padx=8, pady=6)

    status_label.grid(row=5, column=0, columnspan=2, sticky="w", padx=8, pady=(8, 8))
    run_button.grid(row=5, column=2, padx=8, pady=(8, 8))

    def run_clicked() -> None:
        alleles_path = alleles_entry.get().strip()
        genome_dir = genome_entry.get().strip()
        output_path = output_entry.get().strip()

        if not alleles_path:
            messagebox.showerror("Missing input", "Allele FASTA file or directory is required.")
            return
        if not genome_dir:
            messagebox.showerror("Missing input", "Genome directory is required.")
            return
        if not output_path:
            messagebox.showerror("Missing input", "Output CSV path is required.")
            return

        alleles = Path(alleles_path)
        if not alleles.exists():
            messagebox.showerror("Invalid input", f"Allele path not found: {alleles}")
            return
        genomes = Path(genome_dir)
        if not genomes.is_dir():
            messagebox.showerror("Invalid input", f"Genome directory not found: {genomes}")
            return

        try:
            threads = int(threads_entry.get().strip())
        except ValueError:
            messagebox.showerror("Invalid input", "Threads must be an integer.")
            return
        if threads < 1:
            messagebox.showerror("Invalid input", "Threads must be >= 1.")
            return

        progress_var.set(0)
        progress_bar["maximum"] = 1
        run_button.config(state=tk.DISABLED)
        status_var.set("Running...")

        def progress_callback(current: int, total: int, sample_id: str) -> None:
            def update() -> None:
                progress_bar["maximum"] = total
                progress_var.set(current)
                status_var.set(f"Running {current}/{total}: {sample_id}")

            root.after(0, update)

        def worker() -> None:
            try:
                run_pipeline(
                    alleles_path=str(alleles),
                    genome_paths=[],
                    genome_dir=str(genomes),
                    output_path=output_path,
                    min_identity=0.95,
                    strict=False,
                    threads=threads,
                    missing_value="NA",
                    novel_fasta=None,
                    tmp_dir=None,
                    keep_temp=False,
                    allow_locus_only=True,
                    max_genomes=None,
                    progress_callback=progress_callback,
                )
            except Exception as exc:  # pylint: disable=broad-except
                error_message = str(exc)

                def show_error(message: str = error_message) -> None:
                    status_var.set("Failed")
                    run_button.config(state=tk.NORMAL)
                    messagebox.showerror("Run failed", message)

                root.after(0, show_error)
            else:
                def show_success() -> None:
                    status_var.set("Done")
                    run_button.config(state=tk.NORMAL)
                    messagebox.showinfo("Completed", "cgMLST run completed.")

                root.after(0, show_success)

        threading.Thread(target=worker, daemon=True).start()

    run_button.config(command=run_clicked)

    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
