select ST_AsText(geometry) from surface_geometry where parent_id=63;
select ST_AsText(envelope) from cityobject where id=2;
select ST_AsText(geometry) from surface_geometry where parent_id in (select lod2_multi_surface_id from citydb_view.thematic_surface where building_id=2 and objectclass_id=33);