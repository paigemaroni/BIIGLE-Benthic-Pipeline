# Raw Data Archive

`Raw/` is reserved for the original, unmodified BIIGLE downloads and other source material retained for provenance.

For the current pipeline, the BIIGLE ZIP exports are manually concatenated before the formal workflow begins. The resulting canonical starting table is placed at:

```text
Sheets/combined_biigle_annotations.csv
```

The contents of `Raw/` are ignored by Git by default because original exports, imagery and video-derived assets can become large. Keep them locally and back them up using an appropriate institutional, project or archival storage system.

Do not edit raw source files in place. If an input must be corrected, document the correction and create a new canonical prerequisite file under `Sheets/`.

## Worked-example raw provenance

The worked example was manually concatenated from these 25 BIIGLE CSV image-annotation report archives:

```text
32643_csv_image_annotation_report.zip
32572_csv_image_annotation_report.zip
32574_csv_image_annotation_report.zip
32575_csv_image_annotation_report.zip
32581_csv_image_annotation_report.zip
32583_csv_image_annotation_report.zip
32584_csv_image_annotation_report.zip
32585_csv_image_annotation_report.zip
32586_csv_image_annotation_report.zip
32622_csv_image_annotation_report.zip
32625_csv_image_annotation_report.zip
32626_csv_image_annotation_report.zip
32627_csv_image_annotation_report.zip
32629_csv_image_annotation_report.zip
32630_csv_image_annotation_report.zip
32633_csv_image_annotation_report.zip
32634_csv_image_annotation_report.zip
32635_csv_image_annotation_report.zip
32636_csv_image_annotation_report.zip
32637_csv_image_annotation_report.zip
32638_csv_image_annotation_report.zip
32639_csv_image_annotation_report.zip
32640_csv_image_annotation_report.zip
32641_csv_image_annotation_report.zip
32642_csv_image_annotation_report.zip
```

These archives are intentionally ignored by Git by default. The downstream pipeline begins with `Sheets/combined_biigle_annotations.csv`.
