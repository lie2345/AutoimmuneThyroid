"""Correct Chromium's 7.995 pt quantization of nominal 8 pt SVG text.

Only Tf font-size operands associated with effective 7.99--8 pt text change.
All text-showing operands, positioning matrices and non-text operations remain
identical. Raw browser PDFs and an operation-by-operation ledger are retained.
"""
from pathlib import Path
from collections import Counter, defaultdict
import hashlib
import json
import math
import shutil
import sys
import re

import pdfplumber
import pypdfium2 as pdfium
from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf.generic import ContentStream, DecodedStreamObject

TARGET_PT = 8.00005
DPI = 600


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tf_token_spans(data, operations):
    """Lex PDF content while skipping strings/comments, then match Tf operators."""
    whitespace = b'\x00\t\n\x0c\r '
    delimiters = b'()<>[]{}/%'
    tokens=[]
    i=0
    while i < len(data):
        c=data[i]
        if c in whitespace:
            i+=1;continue
        if c==ord('%'):
            while i<len(data) and data[i] not in b'\r\n': i+=1
            continue
        if c==ord('('):
            i+=1;depth=1
            while depth:
                assert i<len(data),'Unterminated PDF string'
                if data[i]==ord('\\'): i+=2;continue
                if data[i]==ord('('): depth+=1
                if data[i]==ord(')'): depth-=1
                i+=1
            continue
        if c==ord('<') and data[i:i+2]!=b'<<':
            end=data.find(b'>',i+1);assert end>=0;i=end+1;continue
        start=i
        if data[i:i+2] in (b'<<',b'>>'):
            i+=2
        elif c==ord('/'):
            i+=1
            while i<len(data) and data[i] not in whitespace+delimiters: i+=1
        elif c in delimiters:
            i+=1
        else:
            while i<len(data) and data[i] not in whitespace+delimiters: i+=1
        tokens.append((data[start:i],start,i))
    # These frozen vector figures contain no inline image syntax.
    assert not any(t[0]==b'BI' for t in tokens)
    lexed=[]
    for i,t in enumerate(tokens):
        if t[0]==b'Tf':
            assert i>=2 and tokens[i-2][0].startswith(b'/')
            number=tokens[i-1]
            assert re.fullmatch(rb'[+-]?(?:\d+\.?\d*|\.\d+)',number[0])
            lexed.append((tokens[i-2][0],number))
    parsed=[(i,args) for i,(args,op) in enumerate(operations) if op==b'Tf']
    assert len(lexed)==len(parsed)
    out={}
    for (name,(token,start,end)),(index,args) in zip(lexed,parsed):
        assert name.decode()==str(args[0]) and float(token)==float(args[1])
        out[index]=(start,end,token)
    return out


def replace_tf_tokens(raw_data, spans, changes):
    replacements=[]
    for change in changes:
        start,end,old=spans[change['operation_index']]
        new=f"{change['Tf_after']:.12f}".rstrip('0').rstrip('.').encode()
        change['raw_token_before']=old.decode()
        change['raw_token_after']=new.decode()
        replacements.append((start,end,old,new))
    revised=raw_data
    for start,end,old,new in reversed(sorted(replacements)):
        assert revised[start:end]==old
        revised=revised[:start]+new+revised[end:]
    # Reverse using new byte offsets; everything outside the Tf values is exact.
    delta=0;reverse=[]
    for start,end,old,new in sorted(replacements):
        reverse.append((start+delta,start+delta+len(new),old,new))
        delta+=len(new)-len(old)
    recovered=revised
    for start,end,old,new in reversed(reverse):
        assert recovered[start:end]==new
        recovered=recovered[:start]+old+recovered[end:]
    assert recovered==raw_data
    return revised


def font_streams(reader):
    result={}
    for name,reference in reader.pages[0]['/Resources']['/Font'].items():
        font=reference.get_object()
        descendant=font.get('/DescendantFonts',[font])[0].get_object()
        descriptor=descendant['/FontDescriptor'].get_object()
        streams={key:hashlib.sha256(descriptor[key].get_object().get_data()).hexdigest()
                 for key in ('/FontFile','/FontFile2','/FontFile3') if key in descriptor}
        assert streams and 'Arial' in font['/BaseFont']
        result[str(name)]={'font':str(font['/BaseFont']),'embedded_stream_sha256':streams}
    return result


def multiply(first, second):
    a,b,c,d,e,f = first
    g,h,i,j,k,l = second
    return [a*g+b*i,a*h+b*j,c*g+d*i,c*h+d*j,e*g+f*i+k,e*h+f*j+l]


def font_scales(operations):
    matrix=[1.,0.,0.,1.,0.,0.]
    stack=[]
    text_matrix=[1.,0.,0.,1.,0.,0.]
    active_font_index=None
    result=defaultdict(list)
    for index,(operands,operator) in enumerate(operations):
        if operator==b'q':
            stack.append(matrix.copy())
        elif operator==b'Q':
            matrix=stack.pop()
        elif operator==b'cm':
            matrix=multiply(list(map(float,operands)),matrix)
        elif operator==b'BT':
            text_matrix=[1.,0.,0.,1.,0.,0.]
        elif operator==b'Tm':
            text_matrix=list(map(float,operands))
        elif operator==b'Tf':
            active_font_index=index
        elif operator in (b'Tj',b'TJ',b"'",b'"'):
            assert active_font_index is not None
            combined=multiply(text_matrix,matrix)
            vertical_scale=math.hypot(combined[2],combined[3])
            result[active_font_index].append(vertical_scale)
    assert not stack
    return result


def normalize_all(root,refresh_existing_outputs=False):
    root=Path(root).resolve()
    validation=root/'validation'
    figure_dir=root/'figures'
    invariant_figure1={str(p.relative_to(root)):sha(p) for p in
        [figure_dir/'Fig1.pdf',figure_dir/'Fig1.png',root/'upload/main_figures/Fig1.tif'] if p.exists()}
    records=[]
    changed_numbers=[]
    for number in range(1,5):
        final=figure_dir/f'Fig{number}.pdf'
        raw=validation/f'main_figure_{number}_raw_browser.pdf'
        if not refresh_existing_outputs or not raw.exists():
            shutil.copy2(final,raw)
        reader=PdfReader(raw)
        before_font_streams=font_streams(reader)
        assert len(reader.pages)==1
        raw_data=reader.pages[0]['/Contents'].get_object().get_data()
        stream=ContentStream(reader.pages[0].get_contents(),reader)
        before_ops=stream.operations
        token_spans=tf_token_spans(raw_data,before_ops)
        scales=font_scales(before_ops)
        with pdfplumber.open(raw) as document:
            before_chars=document.pages[0].chars
            before_sizes=[c['size'] for c in before_chars]
            before_text=[c['text'] for c in before_chars]
            before_fonts=[c['fontname'] for c in before_chars]
            before_box=list(document.pages[0].bbox)
        changes=[]
        for index,scale_list in scales.items():
            operands,operator=before_ops[index]
            assert operator==b'Tf'
            original_size=float(operands[1])
            effective=[original_size*scale for scale in scale_list]
            if min(effective)>=7.99 and max(effective)<8.0:
                corrected_size=TARGET_PT/min(scale_list)
                corrected_effective=[corrected_size*scale for scale in scale_list]
                assert 8.0<=min(corrected_effective)<=max(corrected_effective)<=8.0001
                changes.append({'operation_index':index,'font_resource':str(operands[0]),
                                'Tf_before':original_size,'Tf_after':corrected_size,
                                'effective_before':sorted(set(effective)),
                                'effective_after':sorted(set(corrected_effective)),
                                'effective_vertical_scales':sorted(set(scale_list))})
        if changes:
            assert number in (2,3,4),'Figure 1 must be unaffected'
            writer=PdfWriter()
            writer.clone_document_from_reader(reader)
            revised_data=replace_tf_tokens(raw_data,token_spans,changes)
            updated=DecodedStreamObject()
            updated.set_data(revised_data)
            writer.pages[0].replace_contents(updated)
            temp=validation/f'main_figure_{number}_normalized.tmp.pdf'
            with temp.open('wb') as output:
                writer.write(output)
            reparsed_reader=PdfReader(temp)
            assert reparsed_reader.pages[0]['/Contents'].get_object().get_data()==revised_data
            after_ops=ContentStream(reparsed_reader.pages[0].get_contents(),reparsed_reader).operations
            assert len(after_ops)==len(before_ops)
            changed_indices={c['operation_index'] for c in changes}
            for index,((old_args,old_op),(new_args,new_op)) in enumerate(zip(before_ops,after_ops)):
                assert old_op==new_op
                if index in changed_indices:
                    assert old_op==b'Tf' and old_args[0]==new_args[0] and old_args[1]!=new_args[1]
                else:
                    assert old_args==new_args,(number,index,old_op,old_args,new_args)
            temp.replace(final)
            changed_numbers.append(number)
        with pdfplumber.open(final) as document:
            page=document.pages[0]
            after_chars=page.chars
            assert [c['text'] for c in after_chars]==before_text
            assert [c['fontname'] for c in after_chars]==before_fonts
            assert list(page.bbox)==before_box
            assert len(after_chars)==len(before_chars)
            after_sizes=[c['size'] for c in after_chars]
            assert min(after_sizes)>=8.0 and max(after_sizes)<=12.0,(number,min(after_sizes),max(after_sizes))
        assert before_font_streams==font_streams(PdfReader(final))
        records.append({'number':number,'raw_browser_pdf':str(raw.relative_to(root)),
                        'raw_browser_sha256':sha(raw),'final_pdf_sha256':sha(final),
                        'operation_count':len(before_ops),'changed_Tf_operands':changes,
                        'changed_Tf_count':len(changes),'all_other_operands_and_operators_identical':True,
                        'direct_numeric_operand_comparison':True,
                        'decoded_stream_reverse_substitution_byte_identical':True,
                        'text_strings_fonts_character_count_and_page_box_identical':True,
                        'embedded_Arial_font_streams_identical':before_font_streams,
                        'character_count':len(after_chars),'minimum_pt_before':min(before_sizes),
                        'minimum_pt_after':min(after_sizes),'maximum_pt_after':max(after_sizes),
                        'strict_8_to_12_pass':True})
    if refresh_existing_outputs:
        export_path=validation/'main_figure_export_audit.json'
        export=json.loads(export_path.read_text())
        for number in changed_numbers:
            pdf=figure_dir/f'Fig{number}.pdf'
            doc=pdfium.PdfDocument(pdf)
            rendered=doc[0].render(scale=DPI/72).to_pil().convert('RGB')
            tif=root/'upload/main_figures'/f'Fig{number}.tif'
            rendered.save(tif,format='TIFF',compression='tiff_lzw',dpi=(DPI,DPI))
            rendered.save(figure_dir/f'Fig{number}.png',dpi=(DPI,DPI))
            preview=rendered.copy();preview.thumbnail((1100,1300),Image.Resampling.LANCZOS)
            preview.save(validation/f'main_figure_{number}_preview.png')
            with Image.open(tif) as check:
                assert check.n_frames==1 and check.mode=='RGB' and tuple(check.tag_v2[258])==(8,8,8)
                assert check.tag_v2[259]==5 and 338 not in check.tag_v2
                assert tuple(check.info['dpi'])==(600.,600.) and check.tobytes()==rendered.tobytes()
                assert tif.stat().st_size<10_000_000
            with pdfplumber.open(pdf) as document:
                page=document.pages[0]
                chars=[c for c in page.chars if c['text'].strip()]
                sizes=Counter(round(c['size'],8) for c in chars)
                families=Counter(c['fontname'] for c in chars)
                assert min(sizes)>=8. and max(sizes)<=12.
                assert all('Arial' in name for name in families)
                bounds=[min(c['x0'] for c in chars),min(c['top'] for c in chars),max(c['x1'] for c in chars),max(c['bottom'] for c in chars)]
                assert bounds[0]>=0 and bounds[1]>=0 and bounds[2]<=page.width and bounds[3]<=page.height
            record=export['figures'][number-1]
            record.update({'sha256':sha(tif),'bytes':tif.stat().st_size,'font_sizes_pt':dict(sizes),
                           'actual_pdf_fonts':dict(families),'text_bounds_pt':bounds})
        for record in export['figures']:
            record['checks'].pop('font_size_8_to_12pt_with_0_02pt_renderer_tolerance',None)
            record['checks']['font_size_8_to_12pt']=True
            record['checks']['font_size_8_to_12pt_strict']=True
        export['bounded_renderer_difference']='SVG canvas remains exact. Chromium paper-size quantization is retained; nominal 8 pt text initially at 7.995 pt is corrected to effective 8.00005 pt using only Tf operands. See main_figure_pdf_text_normalization.json.'
        export_path.write_text(json.dumps(export,ensure_ascii=False,indent=2))
    preserved={path:sha(root/path)==digest for path,digest in invariant_figure1.items()}
    assert all(preserved.values())
    audit={'reason':'Correct raw-browser nominal 8 pt text emitted as 7.995 pt; strict journal minimum is 8 pt.',
           'target_effective_pt':TARGET_PT,'figures':records,'changed_figures':changed_numbers,
           'figure1_files_byte_identical':preserved,
           'all_Tm_cm_text_showing_and_nontext_operations_preserved':True,
           'raw_decoded_stream_changed_only_at_validated_Tf_numeric_tokens':True,
           'precision_repair_note':'ContentStream reserialization was replaced with raw token substitution to preserve all original decimal precision; final comparisons use numeric equality and complete reverse-substitution byte identity.',
           'source_SVG_text_sizes_geometry_and_strings_changed':False,
           'character_advance_note':'Within affected text objects, natural glyph advances follow the corrected font size. All explicitly stored text positioning and transformation operands remain identical.'}
    (validation/'main_figure_pdf_text_normalization.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
    for record in records:
        print('Fig',record['number'],'changed Tf',record['changed_Tf_count'],'strict min',record['minimum_pt_after'])
    return audit


if __name__=='__main__':
    normalize_all(Path(__file__).resolve().parents[1],refresh_existing_outputs=True)
