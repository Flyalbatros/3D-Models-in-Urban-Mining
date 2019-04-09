select c.id, ST_AsText(c.envelope) from cityobject c where c.id in 
(select b.id from building b where b.building_parent_id is null) 
union 
select c.id, ST_AsText(c.envelope) from cityobject c where c.objectclass_id = 25 and c.id not in
(select distinct(building_parent_id) from building where building_parent_id is not null);
